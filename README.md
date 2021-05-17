# Suloea.lv2

This is an LV2 plugin wrapper around the excellent [Aeolus] sound engine, by F. Adriaensen.
It only contains a single division and is intended to be used for a single manual.
The original [Aeolus] reverb is included and can be tweaked in part through control ports.
No LV2 gui is done but a generic interface works well for this use case.
A `modgui` is available for use on [MOD devices].

![screenshot][screenshot]

## License

The LV2 plugin wrapper is distributed under the terms of the GPL version 3; see the `LICENSE.md` file.
[Aeolus] is copyright by Fons Adriaensen and Hans Fugal, and distributed under the terms of the GPL version 3; see the information in the [vendored Aeolus] directory.

## Releases

A continuous build is available on the [Releases] page.
Just copy the plugins somewhere your host will find them.

## Compilation

You need to have `cmake` installed, as well as a reasonable compiler.
The procedure is as follows, from the source directory
```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```
followed optionnally by a `sudo make install`.

[Aeolus]: http://kokkinizita.linuxaudio.org/linuxaudio/aeolus/index.html
[vendored Aeolus]: third_party/aeolus-0.9.9
[MOD devices]: https://moddevices.com/
[screenshot]: lv2/modgui/screenshot-suloea.png