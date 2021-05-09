/*
  Toccata LV2 plugin

  Copyright 2020, Paul Ferrand <paul@ferrand.cc>

  This file was based on skeleton and example code from the LV2 plugin
  distribution available at http://lv2plug.in/

  The LV2 sample plugins have the following copyright and notice, which are
  extended to the current work:
  Copyright 2011-2016 David Robillard <d@drobilla.net>
  Copyright 2011 Gabriel M. Beddingfield <gabriel@teuton.org>
  Copyright 2011 James Morris <jwm.art.net@gmail.com>

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THIS SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "lv2/atom/forge.h"
#include "lv2/atom/util.h"
#include "lv2/buf-size/buf-size.h"
#include "lv2/core/lv2.h"
#include "lv2/core/lv2_util.h"
#include "lv2/midi/midi.h"
#include "lv2/options/options.h"
#include "lv2/parameters/parameters.h"
#include "lv2/patch/patch.h"
#include "lv2/state/state.h"
#include "lv2/urid/urid.h"
#include "lv2/worker/worker.h"
#include "lv2/time/time.h"
#include "lv2/log/logger.h"
#include "lv2/log/log.h"

#include "division.h"
#include "asection.h"
#include "scales.h"
#include "reverb.h"
#include "stops.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <memory>
#include <vector>

#define AEOLUS_PERIOD 64
#define SULOEA_URI "https://github.com/paulfd/suloea"
#define CHANNEL_MASK 0x0F
#define NOTE_ON 0x90
#define NOTE_OFF 0x80
#define MIDI_CHANNEL(byte) (byte & CHANNEL_MASK)
#define MIDI_STATUS(byte) (byte & ~CHANNEL_MASK)
#define MAX_BLOCK_SIZE 8192
#define MAX_PATH_SIZE 1024
#define NUM_VOICES 256
#define UNUSED(x) (void)(x)

enum {
    INPUT_PORT = 0,
    LEFT_BUFFER,
    RIGHT_BUFFER,
    FREEWHEEL_PORT
};

struct SuloeaPlugin
{
    // Features
    LV2_URID_Map* map;
    LV2_URID_Unmap* unmap;
    LV2_Log_Log* log;

    // Ports
    const LV2_Atom_Sequence* input_port;
    float *output_buffers[2];
    const float *freewheel_port;

    // Atom forge
    LV2_Atom_Forge forge; ///< Forge for writing atoms in run thread
    LV2_Atom_Forge_Frame notify_frame; ///< Cached for worker replies

    // Logger
    LV2_Log_Logger logger;

    // URIs
    LV2_URID midi_event_uri;
    LV2_URID options_interface_uri;
    LV2_URID max_block_length_uri;
    LV2_URID nominal_block_length_uri;
    LV2_URID sample_rate_uri;
    LV2_URID atom_object_uri;
    LV2_URID atom_float_uri;
    LV2_URID atom_int_uri;
    LV2_URID atom_urid_uri;
    LV2_URID atom_string_uri;
    LV2_URID atom_bool_uri;

    std::vector<Addsynth> synths;
    Reverb reverb;
    std::unique_ptr<Asection> asection;
    std::unique_ptr<Division> division;
    unsigned char keymap[NNOTES];
    float W[PERIOD];
    float X[PERIOD];
    float Y[PERIOD];
    float Z[PERIOD];
    float R[PERIOD];

    bool activated;
    int max_block_size;
    double sample_rate;

    static std::unique_ptr<SuloeaPlugin> instantiate(double rate, 
        const char* path, const LV2_Feature* const* features);
    // Logic taken from Aeolus
    void key_off (int n, int b)
    {
        keymap [n] &= ~b;
        keymap [n] |= 128;
    }

    // Logic taken from Aeolus
    void key_on (int n, int b)
    {
        keymap [n] |= b | 128;
    }

    // Logic taken from Aeolus (proc_keys1() and proc_keys2() combined)
    void proc_keys()
    {
        int m, n;

        for (n = 0; n < NNOTES; n++)
        {
            m = keymap[n];
            if (m & 128)
            {
                m &= 127;
                keymap[n] = m;
                division->update(n, m);
            }
        }

        division->update(keymap);
    }

    void process_midi_event(const LV2_Atom_Event* ev);
    void check_freewheeling();
    void map_required_uris();
    void process_output(uint32_t sample_count);
};

void SuloeaPlugin::map_required_uris()
{
    midi_event_uri = map->map(map->handle, LV2_MIDI__MidiEvent);
    max_block_length_uri = map->map(map->handle, LV2_BUF_SIZE__maxBlockLength);
    nominal_block_length_uri = map->map(map->handle, LV2_BUF_SIZE__nominalBlockLength);
    sample_rate_uri = map->map(map->handle, LV2_PARAMETERS__sampleRate);
    atom_float_uri = map->map(map->handle, LV2_ATOM__Float);
    atom_int_uri = map->map(map->handle, LV2_ATOM__Int);
    atom_bool_uri = map->map(map->handle, LV2_ATOM__Bool);
    atom_string_uri = map->map(map->handle, LV2_ATOM__String);
    atom_urid_uri = map->map(map->handle, LV2_ATOM__URID);
}

std::unique_ptr<SuloeaPlugin> SuloeaPlugin::instantiate(double rate, 
        const char* path, const LV2_Feature* const* features)
{
    LV2_Options_Option* options = NULL;
    bool supports_bounded_block_size = false;
    bool options_has_block_size = false;
    bool supports_fixed_block_size = false;

    // Allocate and initialise instance structure.
    std::unique_ptr<SuloeaPlugin> self { new SuloeaPlugin };
    if (!self)
        return {};

    // Set defaults
    self->max_block_size = MAX_BLOCK_SIZE;
    self->sample_rate = rate;
    self->activated = false;

    // Get the features from the host and populate the structure
    for (const LV2_Feature* const* f = features; *f; f++) {
        void *data = (**f).data;

        if (!strcmp((**f).URI, LV2_URID__map))
            self->map = (LV2_URID_Map *)data;

        if (!strcmp((**f).URI, LV2_URID__unmap))
            self->unmap = (LV2_URID_Unmap *)data;

        if (!strcmp((**f).URI, LV2_BUF_SIZE__boundedBlockLength))
            supports_bounded_block_size = true;

        if (!strcmp((**f).URI, LV2_BUF_SIZE__fixedBlockLength))
            supports_fixed_block_size = true;

        if (!strcmp((**f).URI, LV2_OPTIONS__options))
            options = (LV2_Options_Option *)data;

        if (!strcmp((**f).URI, LV2_LOG__log))
            self->log = (LV2_Log_Log *)data;
    }

    // Setup the loggers
    lv2_log_logger_init(&self->logger, self->map, self->log);

    // The map feature is required
    if (!self->map) {
        lv2_log_error(&self->logger, "Map feature not found, aborting...\n");
        return {};
    }

    // Map the URIs we will need
    self->map_required_uris();

    // Initialize the forge
    lv2_atom_forge_init(&self->forge, self->map);

    // Check the options for the block size and sample rate parameters
    if (options) {
        for (const LV2_Options_Option* opt = options; opt->value || opt->key; ++opt) {
            if (opt->key == self->sample_rate_uri) {
                if (opt->type != self->atom_float_uri) {
                    lv2_log_warning(&self->logger, "Got a sample rate but the type was wrong\n");
                    continue;
                }
                self->sample_rate = *(float*)opt->value;
            } else if (opt->key == self->max_block_length_uri) {
                if (opt->type != self->atom_int_uri) {
                    lv2_log_warning(&self->logger, "Got a max block size but the type was wrong\n");
                    continue;
                }
                self->max_block_size = *(int*)opt->value;
                options_has_block_size = true;
            } else if (opt->key == self->nominal_block_length_uri) {
                if (opt->type != self->atom_int_uri) {
                    lv2_log_warning(&self->logger, "Got a nominal block size but the type was wrong\n");
                    continue;
                }
                self->max_block_size = *(int*)opt->value;
                options_has_block_size = true;
            }
        }
    } else {
        lv2_log_warning(&self->logger,
            "No option array was given upon instantiation; will use default values\n.");
    }

    // We need _some_ information on the block size
    if (!supports_bounded_block_size && !supports_fixed_block_size && !options_has_block_size) {
        lv2_log_error(&self->logger,
            "Bounded block size not supported and options gave no block size, aborting...\n");
        return {};
    }

    memset(self->keymap, 0, NNOTES * sizeof(unsigned char));
    self->synths.resize(stopList.size());
    std::string stopPath { path };
    stopPath += "/stops";
    for (unsigned i = 0; i < stopList.size(); ++i) {
        const StopDescription& stop = stopList[i];
        Addsynth& synth = self->synths[i];
        strcpy(synth._filename, stop.filename);
        synth._pan = stop.pan;
        synth._del = stop.del;
        synth.load(stopPath.c_str());
    }

    self->asection.reset(new Asection(self->sample_rate));
    self->asection->set_size(0.075f); // In init_audio (audio.cc)
    self->division.reset(new Division(self->asection.get(), self->sample_rate));
    self->reverb.init(self->sample_rate);

    // Use min/max instead of the original range. Why not?
    for (unsigned i = 0; i < stopList.size(); ++i) {
        const StopDescription& stop = stopList[i];
        Addsynth& synth = self->synths[i];
        std::unique_ptr<Rankwave> wave { new Rankwave(NOTE_MIN, NOTE_MAX) };
        wave->gen_waves(&synth, self->sample_rate, 440.0f, scales[4]._data);
        self->division->set_rank(i, wave.release(), stop.pan, stop.del);
        
        // Enable the rank
        self->division->set_rank_mask(i, 128);
    }

    self->division->set_div_mask(1); // This is also the value that gets entered in key on/off
                                     // Probably used for coupling logic overall in the original
                                     // Aeolus
    return self;
}


static void
connect_port(LV2_Handle instance,
    uint32_t port,
    void* data)
{
    SuloeaPlugin* self = (SuloeaPlugin*)instance;
    switch (port) {
    case INPUT_PORT:
        self->input_port = (const LV2_Atom_Sequence*)data;
        break;
    case LEFT_BUFFER:
        self->output_buffers[0] = (float *)data;
        break;
    case RIGHT_BUFFER:
        self->output_buffers[1] = (float *)data;
        break;
    case FREEWHEEL_PORT:
        self->freewheel_port = (const float *)data;
        break;
    default:
        break;
    }
}

static LV2_Handle
instantiate(const LV2_Descriptor* descriptor,
    double rate,
    const char* path,
    const LV2_Feature* const* features)
{
    UNUSED(descriptor);
    auto self = SuloeaPlugin::instantiate(rate, path, features);
    return (LV2_Handle)self.release();
}

static void
cleanup(LV2_Handle instance)
{
    SuloeaPlugin* self = (SuloeaPlugin*)instance;
    free(self);
}

static void
activate(LV2_Handle instance)
{
    SuloeaPlugin* self = (SuloeaPlugin*)instance;
    self->activated = true;
}

static void
deactivate(LV2_Handle instance)
{
    SuloeaPlugin* self = (SuloeaPlugin*)instance;
    self->activated = false;
}

void SuloeaPlugin::process_midi_event(const LV2_Atom_Event* ev)
{
    const uint8_t* const msg = (const uint8_t*)(ev + 1);
    switch (lv2_midi_message_type(msg)) {
    case LV2_MIDI_MSG_NOTE_ON:
        if (msg[2] == 0)
            goto noteoff; // 0 velocity note-ons should be forbidden but just in case...
        key_on((int)msg[1] - NOTE_MIN, 1);
        break;
    case LV2_MIDI_MSG_NOTE_OFF: noteoff:
        key_off((int)msg[1] - NOTE_MIN, 1);
        break;
    default:
        break;
    }
}

void SuloeaPlugin::check_freewheeling()
{
    if (*freewheel_port > 0)
    {
        // NOOP
    }
    else
    {
        // NOOP
    }
}

void SuloeaPlugin::process_output(uint32_t sample_count)
{
    // Remember to set reverb parameters (see proc_synth in audio.cc)

    float* out [2] { output_buffers[0], output_buffers[1] };
    for (uint32_t k = 0; k < sample_count; k += PERIOD)
    {
        memset(W, 0, PERIOD * sizeof(float));
        memset(X, 0, PERIOD * sizeof(float));
        memset(Y, 0, PERIOD * sizeof(float));
        memset(Z, 0, PERIOD * sizeof(float));
        memset(R, 0, PERIOD * sizeof(float));

        division->process();
        // Volume is a default value in audio.cc init_audio() _audiopar[VOLUME]
        asection->process(0.32f, W, X, Y, R);
        reverb.process(PERIOD, 0.32f, R, W, X, Y, Z);
        // Check proc_synth in audio.cc for the ambisonic version

        for (unsigned j = 0; j < PERIOD; j++)
        {
            // Default value in audio.cc init_audio() _audiopar[STPOSIT], 
            // seems like a stereo position
            out[0][j] = 
                W[j] + 0.5f * X[j] + Y[j];
            out[1][j] =
                W[j] + 0.5f * X[j] - Y[j];
        }

        out[0] += PERIOD;
        out[1] += PERIOD;
    } 
}

static void
run(LV2_Handle instance, uint32_t sample_count)
{
    UNUSED(sample_count);
    SuloeaPlugin* self = (SuloeaPlugin*)instance;
    if (!self->activated)
        return;

    if (!self->input_port)
        return;

    LV2_ATOM_SEQUENCE_FOREACH(self->input_port, ev)
    {
        // If the received atom is an object/patch message
        if (ev->body.type == self->atom_object_uri) {
            const LV2_Atom_Object* obj = (const LV2_Atom_Object*)&ev->body;
            lv2_log_warning(&self->logger, "Got an Object atom but it was not supported.\n");
            if (self->unmap)
                lv2_log_warning(&self->logger,
                    "Object URI: %s\n",
                    self->unmap->unmap(self->unmap->handle, obj->body.otype));
            continue;
            // Got an atom that is a MIDI event
        } else if (ev->body.type == self->midi_event_uri) {
            self->process_midi_event(ev);
        }
    }

    self->proc_keys();
    self->check_freewheeling();
    self->process_output(sample_count);
}

static uint32_t
lv2_get_options(LV2_Handle instance, LV2_Options_Option* options)
{
    UNUSED(instance);
    UNUSED(options);
    // We have no options
    return LV2_OPTIONS_ERR_UNKNOWN;
}

static uint32_t
lv2_set_options(LV2_Handle instance, const LV2_Options_Option* options)
{
    SuloeaPlugin* self = (SuloeaPlugin*)instance;

    // Update the block size and sample rate as needed
    for (const LV2_Options_Option* opt = options; opt->value || opt->key; ++opt) {
        if (opt->key == self->sample_rate_uri) {
            if (opt->type != self->atom_float_uri) {
                lv2_log_warning(&self->logger, "Got a sample rate but the type was wrong\n");
                continue;
            }
            self->sample_rate = *(float*)opt->value;
            // sfizz_set_sample_rate(self->synth, self->sample_rate);
        } else if (opt->key == self->nominal_block_length_uri) {
            if (opt->type != self->atom_int_uri) {
                lv2_log_warning(&self->logger, "Got a nominal block size but the type was wrong\n");
                continue;
            }
            self->max_block_size = *(int*)opt->value;
        }
    }
    return LV2_OPTIONS_SUCCESS;
}

static const void*
extension_data(const char* uri)
{
    static const LV2_Options_Interface options = { lv2_get_options, lv2_set_options };
    // Advertise the extensions we support
    if (!strcmp(uri, LV2_OPTIONS__interface))
        return &options;

    return NULL;
}

static const LV2_Descriptor descriptor = {
    SULOEA_URI,
    instantiate,
    connect_port,
    activate,
    run,
    deactivate,
    cleanup,
    extension_data
};

LV2_SYMBOL_EXPORT
const LV2_Descriptor*
lv2_descriptor(uint32_t index)
{
    switch (index) {
    case 0:
        return &descriptor;
    default:
        return NULL;
    }
}