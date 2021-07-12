/*
  SPDX-License-Identifier: GPL-3.0-only
  
  Suolea LV2 plugin wrapper

  Copyright 2021, Paul Ferrand <paul@ferrand.cc>

  This file is licensed under the GNU General Public License version 3.
  For a full text of the license, see LICENSE.md.

  This file was based on skeleton and example code from the LV2 plugin
  distribution available at http://lv2plug.in/

  The LV2 sample plugins have the following copyright and notice:
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

#include <cmath>
#include <string>
#include <memory>
#include <vector>

#define SULOEA_URI "https://github.com/paulfd/suloea"
#define SULOEA__retune SULOEA_URI ":" "retune"
#define MAX_BLOCK_SIZE 8192
#define UNUSED(x) (void)(x)
#define ASECTION_REVERB 4 // see asection.h

namespace Default {
    constexpr float volume { -10.0f };
    constexpr float reverb_time { 4.0f };
    constexpr float reverb_delay { 50.0f };
    constexpr float reverb_amount { -10.0f };
    constexpr float frequency { 440.0f };
    constexpr int scale_index { 4 };
}

enum {
    INPUT_PORT = 0,
    LEFT_BUFFER,
    RIGHT_BUFFER,
    VOLUME_PORT,
    REVERB_PORT,
    DELAY_PORT,
    TIME_PORT,
    REVERB_AMOUNT_PORT,
    POSITION_PORT,
    FREQUENCY_PORT,
    SCALE_PORT,
    RETUNING_PORT,
    PORT_SENTINEL
};

struct SuloeaPlugin
{
    // Features
    LV2_URID_Map* map;
    LV2_URID_Unmap* unmap;
    LV2_Log_Log* log;
    LV2_Worker_Schedule *worker {};

    // Ports
    const LV2_Atom_Sequence* input_port;
    float *output_buffers[2];
    const float *volume_port;
    const float *reverb_port;
    const float *delay_port;
    const float *time_port;
    const float *reverb_amount_port;
    const float *position_port;
    const float *scale_port;
    const float *frequency_port;
    float *retuning_port;
    std::vector<const float*> stop_ports;

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
    LV2_URID atom_double_uri;
    LV2_URID atom_int_uri;
    LV2_URID atom_urid_uri;
    LV2_URID atom_string_uri;
    LV2_URID atom_bool_uri;
    LV2_URID retune_uri;

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
    unsigned read_idx { 0 };

    bool activated;
    bool retuning { false };
    int max_block_size;
    double sample_rate;
    bool reverb_active { true };
    float reverb_time { Default::reverb_time };
    float reverb_delay { Default::reverb_delay };
    float reverb_amount { Default::reverb_amount };
    float volume { Default::volume };
    int scale_index { Default::scale_index };
    float frequency { Default::frequency };
    float gain;
    float reverb_amount_gain;
    std::vector<bool> active_stops;

    static std::unique_ptr<SuloeaPlugin> instantiate(double rate, 
        const char* path, const LV2_Feature* const* features);
    void key_off (int n, int b);
    void key_on (int n, int b);
    void proc_keys();
    void retune_stops();
    void init();
    void update_stops();
    void update_parameters();
    void process_midi_event(const LV2_Atom_Event* ev);
    void map_required_uris();
    void clear_internal_buffers();
    void process_one_period();
    void update_gain();
    void update_reverb_delay();
    void update_reverb_time();
    void update_reverb_amount();
};

void SuloeaPlugin::map_required_uris()
{
    midi_event_uri = map->map(map->handle, LV2_MIDI__MidiEvent);
    max_block_length_uri = map->map(map->handle, LV2_BUF_SIZE__maxBlockLength);
    nominal_block_length_uri = map->map(map->handle, LV2_BUF_SIZE__nominalBlockLength);
    sample_rate_uri = map->map(map->handle, LV2_PARAMETERS__sampleRate);
    atom_float_uri = map->map(map->handle, LV2_ATOM__Float);
    atom_double_uri = map->map(map->handle, LV2_ATOM__Double);
    atom_int_uri = map->map(map->handle, LV2_ATOM__Int);
    atom_bool_uri = map->map(map->handle, LV2_ATOM__Bool);
    atom_string_uri = map->map(map->handle, LV2_ATOM__String);
    atom_urid_uri = map->map(map->handle, LV2_ATOM__URID);
    retune_uri = map->map(map->handle, SULOEA__retune);
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
        const char *uri = (**f).URI;

        if (!strcmp(uri, LV2_URID__map))
            self->map = (LV2_URID_Map *)data;

        if (!strcmp(uri, LV2_URID__unmap))
            self->unmap = (LV2_URID_Unmap *)data;

        if (!strcmp(uri, LV2_BUF_SIZE__boundedBlockLength))
            supports_bounded_block_size = true;

        if (!strcmp(uri, LV2_BUF_SIZE__fixedBlockLength))
            supports_fixed_block_size = true;

        if (!strcmp(uri, LV2_OPTIONS__options))
            options = (LV2_Options_Option *)data;

        if (!strcmp(uri, LV2_LOG__log))
            self->log = (LV2_Log_Log *)data;

        if (!strcmp(uri, LV2_WORKER__schedule))
            self->worker = (LV2_Worker_Schedule *)data;
    }

    // Setup the loggers
    lv2_log_logger_init(&self->logger, self->map, self->log);

    // The map feature is required
    if (!self->map || !self->unmap) {
        lv2_log_error(&self->logger, "Map/unmap feature not found, aborting...\n");
        return {};
    }

    // The worker feature is required
    if (!self->worker)
    {
        lv2_log_error(&self->logger, "Worker feature not found, aborting..\n");
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
                if (opt->type == self->atom_float_uri) {
                    self->sample_rate = *(float*)opt->value;
                } else if (opt->type == self->atom_double_uri) {
                    self->sample_rate = *(double*)opt->value;
                } else if (opt->type == self->atom_int_uri) {
                    // TODO : MOD host sends 1024 as int here?
                    // self->sample_rate = (double)*(int*)opt->value;
                    // lv2_log_warning(&self->logger, "Sample rate: %f\n", self->sample_rate);
                } else {
                    lv2_log_warning(&self->logger,
                          "Got an unhandled sample rate type: %s\n",
                          self->unmap->unmap(self->unmap->handle, opt->type));
                }
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
    for (unsigned i = 0; i < stopList.size(); ++i) {
        self->stop_ports.push_back(nullptr);
        self->active_stops.push_back(false);
    }

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
    self->division.reset(new Division(self->asection.get(), self->sample_rate));
    self->division->set_div_mask(1); // This is also the value that gets entered in key on/off
                                     // Probably used for coupling logic overall in the original
                                     // Aeolus
    self->reverb.init(self->sample_rate);
    self->init();
    return self;
}

void SuloeaPlugin::retune_stops()
{
    for (unsigned i = 0; i < stopList.size(); ++i) {
        const StopDescription& stop = stopList[i];
        Addsynth& synth = synths[i];
        std::unique_ptr<Rankwave> wave { new Rankwave(NOTE_MIN, NOTE_MAX) };
        wave->gen_waves(&synth, sample_rate, frequency, scales[scale_index]._data);
        division->set_rank(i, wave.release(), stop.pan, stop.del);
    }
}

// Logic taken from Aeolus
void SuloeaPlugin::key_off (int n, int b)
{
    keymap [n] &= ~b;
    keymap [n] |= 128;
}

// Logic taken from Aeolus
void SuloeaPlugin::key_on (int n, int b)
{
    keymap [n] |= b | 128;
}

// Logic taken from Aeolus (proc_keys1())
// proc_keys2() is in update_stops()
void SuloeaPlugin::proc_keys()
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
}

void SuloeaPlugin::update_gain()
{
    gain = std::pow(10.0f, volume / 20.0f);
}

void SuloeaPlugin::update_reverb_delay()
{ 
    float delay_in_seconds = reverb_delay * 1e-3f;
    reverb.set_delay(delay_in_seconds);
    asection->set_size(delay_in_seconds);
}

void SuloeaPlugin::update_reverb_time()
{
    reverb.set_t60mf (reverb_time);
    reverb.set_t60lo (reverb_time * 1.50f, 250.0f);
    reverb.set_t60hi (reverb_time * 0.50f, 3e3f);
}

void SuloeaPlugin::init()
{
    clear_internal_buffers();
    update_gain();
    update_reverb_amount();
    update_reverb_delay();
    update_reverb_time();
    retune_stops();

    for (unsigned i = 0; i < active_stops.size(); ++i) {
        if (active_stops[i])
            division->set_rank_mask(i, 128); 
        else
            division->clr_rank_mask(i, 128); 
    }
}

void SuloeaPlugin::update_stops()
{
    for (unsigned i = 0; i < active_stops.size(); ++i) {
        if (stop_ports[i] == nullptr)
            continue;

        bool stop_status = *stop_ports[i] > 0.0f;
        if (stop_status == active_stops[i])
            continue;

        active_stops[i] = stop_status;
        if (stop_status)
            division->set_rank_mask(i, 128); 
        else
            division->clr_rank_mask(i, 128); 
    }

    // Update playing notes
    division->update(keymap);
}

void SuloeaPlugin::update_reverb_amount()
{
    Fparm* reverb_parameters = asection->get_apar();
    if (reverb_active) {
        reverb_parameters[ASECTION_REVERB]._val = std::pow(10.0f, reverb_amount / 10.0f);
    } else {
        reverb_parameters[ASECTION_REVERB]._val = 0.0f;
    }
}

void SuloeaPlugin::update_parameters()
{
    *retuning_port = (float)retuning;

    // From proc_synth() in audio.cc
    if (std::fabs(reverb_delay - *delay_port) > 1) {
        reverb_delay = *delay_port;
        update_reverb_delay();
    }

    if (std::fabs(reverb_time - *time_port) > 0.1f) {
        reverb_time = *time_port;
        update_reverb_time();
    }

    if (std::fabs(volume - *volume_port) > 0.1f) {
        volume = *volume_port;
        update_gain();
    }

    if (reverb_active != (bool)*reverb_port) {
        reverb_active = (bool)*reverb_port;
        update_reverb_amount();
    }

    if (std::fabs(reverb_amount - *reverb_amount_port) > 0.1f) {
        reverb_amount = *reverb_amount_port;
        update_reverb_amount();
    }

    if (scale_index != (int)*scale_port && !retuning) {
        retuning = true;
        scale_index = (int)*scale_port;
        LV2_Atom atom;
        atom.size = 0;
        atom.type = retune_uri;
        if (!(worker->schedule_work(worker->handle,
                                         lv2_atom_total_size((LV2_Atom *)&atom),
                                         &atom) == LV2_WORKER_SUCCESS))
        {
            lv2_log_error(&logger, "[suleoa] There was an issue letting the background worker retune\n");
        }
    }

    if (!retuning && std::fabs(frequency - *frequency_port) > 0.05f ) {
        retuning = true;
        frequency = *frequency_port;
        LV2_Atom atom;
        atom.size = 0;
        atom.type = retune_uri;
        if (!(worker->schedule_work(worker->handle,
                                         lv2_atom_total_size((LV2_Atom *)&atom),
                                         &atom) == LV2_WORKER_SUCCESS))
        {
            lv2_log_error(&logger, "[suleoa] There was an issue letting the background worker retune\n");
        }
    }
}

void SuloeaPlugin::process_midi_event(const LV2_Atom_Event* ev)
{
    const uint8_t* const msg = (const uint8_t*)(ev + 1);
    switch (lv2_midi_message_type(msg)) {
    case LV2_MIDI_MSG_NOTE_ON:
        if (msg[2] == 0)
            goto noteoff; // 0 velocity note-ons should be forbidden but just in case...
        if (msg[1] >= NOTE_MIN && msg[1] <= NOTE_MAX)
            key_on((int)msg[1] - NOTE_MIN, 1);
        break;
    case LV2_MIDI_MSG_NOTE_OFF: noteoff:
        if (msg[1] >= NOTE_MIN && msg[1] <= NOTE_MAX)
            key_off((int)msg[1] - NOTE_MIN, 1);
        break;
    default:
        break;
    }
}

void SuloeaPlugin::clear_internal_buffers()
{
    memset(W, 0, PERIOD * sizeof(float));
    memset(X, 0, PERIOD * sizeof(float));
    memset(Y, 0, PERIOD * sizeof(float));
    memset(Z, 0, PERIOD * sizeof(float));
    memset(R, 0, PERIOD * sizeof(float));
}

void SuloeaPlugin::process_one_period()
{
    clear_internal_buffers();

    if (!retuning)
        division->process();

    // Volume is a default value in audio.cc init_audio() _audiopar[VOLUME]
    asection->process(gain, W, X, Y, R);
    reverb.process(PERIOD, gain, R, W, X, Y, Z);
    // Check proc_synth in audio.cc for the ambisonic version
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
    case VOLUME_PORT:
        self->volume_port = (const float *)data;
        break;
    case REVERB_PORT:
        self->reverb_port = (const float *)data;
        break;
    case DELAY_PORT:
        self->delay_port = (const float *)data;
        break;
    case TIME_PORT:
        self->time_port = (const float *)data;
        break;
    case REVERB_AMOUNT_PORT:
        self->reverb_amount_port = (const float *)data;
        break;
    case POSITION_PORT:
        self->position_port = (const float *)data;
        break;
    case SCALE_PORT:
        self->scale_port = (const float *)data;
        break;
    case FREQUENCY_PORT:
        self->frequency_port = (const float *)data;
        break;
    case RETUNING_PORT:
        self->retuning_port = (float *)data;
        break;
    default: // Afterwards it's a stop
        {
            uint32_t stop_index = port - PORT_SENTINEL;
            if (stop_index < self->stop_ports.size())
                self->stop_ports[stop_index] = (const float*)data;
            else
                lv2_log_error(&self->logger, "Port index out of range: %d\n", port);
        }
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

static void
run(LV2_Handle instance, uint32_t sample_count)
{
    UNUSED(sample_count);
    SuloeaPlugin* self = (SuloeaPlugin*)instance;
    if (!self->activated)
        return;

    if (!self->input_port)
        return;

    const LV2_Atom_Sequence* seq = self->input_port;

    self->update_parameters();
    self->update_stops();

    LV2_Atom_Event* ev = lv2_atom_sequence_begin(&seq->body);
    float* out [2] { self->output_buffers[0], self->output_buffers[1] };
    int64_t k = 0;

    if (self->retuning) {
        // Process midi event blanks
        while (!lv2_atom_sequence_is_end(&seq->body, seq->atom.size, ev)) {
            if (ev->body.type == self->midi_event_uri)
                self->process_midi_event(ev);

            ev = lv2_atom_sequence_next(ev);
        }

        // Clear buffers
        memset(self->output_buffers[0], 0, sample_count * sizeof(float));
        memset(self->output_buffers[1], 0, sample_count * sizeof(float));

        return;
    }

    while (k < sample_count)
    {
        int64_t left_to_read = PERIOD - self->read_idx;
        int64_t left_in_block = (int64_t)sample_count - k;
        int64_t chunk_size = left_to_read == 0 ? 
            std::min((int64_t) PERIOD, left_in_block) : std::min(left_to_read, left_in_block);
        int64_t sentinel = k + chunk_size;

        while (!lv2_atom_sequence_is_end(&seq->body, seq->atom.size, ev) 
            && ev->time.frames < sentinel) {
            if (ev->body.type == self->midi_event_uri)
                self->process_midi_event(ev);

            ev = lv2_atom_sequence_next(ev);
        }

        if (left_to_read == 0) {
            self->proc_keys();
            self->process_one_period();
            self->read_idx = 0;
        }
        
        unsigned& j = self->read_idx;
        for (unsigned i = 0; i < chunk_size; ++j, ++i)
        {
            *out[0]++ = self->W[j] + *self->position_port * self->X[j] + self->Y[j];
            *out[1]++ = self->W[j] + *self->position_port * self->X[j] - self->Y[j];
        }

        k += chunk_size;
    }
}

// This runs in a lower priority thread
static LV2_Worker_Status
work(LV2_Handle instance,
     LV2_Worker_Respond_Function respond,
     LV2_Worker_Respond_Handle handle,
     uint32_t size,
     const void *data)
{
    UNUSED(size);
    SuloeaPlugin *self = (SuloeaPlugin *)instance;
    if (!data) {
        lv2_log_error(&self->logger, "[suloea] Ignoring empty data in the worker thread\n");
        return LV2_WORKER_ERR_UNKNOWN;
    }

    const LV2_Atom *atom = (const LV2_Atom *)data;
    if (atom->type == self->retune_uri) {
        self->retune_stops();
        respond(handle, size, data);
    } else {
        lv2_log_error(&self->logger, "[sfizz] Got an unknown atom in work\n");
        if (self->unmap)
            lv2_log_error(&self->logger,
                          "URI: %s\n",
                          self->unmap->unmap(self->unmap->handle, atom->type));
        return LV2_WORKER_ERR_UNKNOWN;
    }
    return LV2_WORKER_SUCCESS;
}

// This runs in the audio thread
static LV2_Worker_Status
work_response(LV2_Handle instance,
              uint32_t size,
              const void *data)
{
    UNUSED(size);
    SuloeaPlugin *self = (SuloeaPlugin *)instance;

    if (!data)
        return LV2_WORKER_ERR_UNKNOWN;

    const LV2_Atom *atom = (const LV2_Atom *)data;
    if (atom->type == self->retune_uri) {
        self->retuning = false; 
    } else {
        lv2_log_error(&self->logger, "[sfizz] Got an unexpected atom in work response\n");
        if (self->unmap)
            lv2_log_error(&self->logger,
                          "URI: %s\n",
                          self->unmap->unmap(self->unmap->handle, atom->type));
        return LV2_WORKER_ERR_UNKNOWN;
    }

    return LV2_WORKER_SUCCESS;
}

static const void*
extension_data(const char* uri)
{
    static const LV2_Worker_Interface worker = { work, work_response, NULL }; 

    // Advertise the extensions we support
    if (!strcmp(uri, LV2_WORKER__interface))
        return &worker;

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