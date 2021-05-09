#pragma once
#include <vector>
#include "addsynth.h"

struct StopDescription
{
    const char* filename;
    int del;
    char pan;
};

const std::vector<StopDescription> stopList {
    { "bourdon16.ae0", 20, 'C' },
    { "flute8.ae0", 5, 'C' },
    { "I_principal_8.ae0", 10, 'C' },
    { "flute4.ae0", 8, 'R' },
    { "I_principal_4.ae0", 14, 'R' },
    { "flute2.ae0", 10, 'R' },
    { "I_octave_2.ae0", 16, 'L' },
    { "I_octave_1.ae0", 19, 'R' },
    { "I_quinte_513.ae0", 21, 'C' },
    // { "I_quinte_223.ae0", 19, 'R' },
    { "tibia8.ae0", 22, 'L' },
    { "celesta.ae0", 13, 'C' },
    { "I_cymbel.ae0", 4, 'C' },
    { "I_mixtur5fach.ae0", 18, 'L' },
    { "I_trumpet.ae0", 0, 'W' },
};