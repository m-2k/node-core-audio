# Fork Information

The purpose of this fork is to strip away unused npm dependencies and 
leave just the core PortAudio wrapper, for use with my other project
[SKQW](https://github.com/michaelbromley/skqw).

I added a method to the API:

```javascript
var engine = coreAudio.createNewAudioEngine();

// set the frequency (in hertz) by which the engine runs
// the processAudio method. Formerly this was 1000Hz, but
// now it defaults to 60Hz.
engine.setSampleRate( 100 );
```


Node Core Audio
==================

![alt tag](https://nodei.co/npm-dl/node-core-audio.png)

A C++ extension for node.js that gives javascript access to audio buffers and basic audio processing functionality

Right now, it's basically a node.js binding for PortAudio.

NOTE: Looking for help maintaining this repository!

Active contributors:

- [rmedaer](https://github.com/rmedaer)

Installation
=====

```
npm install node-core-audio
```

Basic Usage
=====

Below is the most basic use of the audio engine. We create a new instance of
node-core-audio, and then give it our processing function. The audio engine
will call the audio callback whenever it needs an output buffer to send to
the sound card.

```javascript
// Create a new instance of node-core-audio
const coreAudio = require("node-core-audio")

// Create a new audio engine
var engine = coreAudio.createNewAudioEngine()

/**
 * A processing function that can process / manipulate the incoming audio samples
 * before returning them to the soundcard's output.
 * @param {Array} samples The samples from the input stream.
 * @return {Array} The samples for the output stream.
 */
function processAudio( inputBuffer ) {	
	for (var channel = 0; channel < inputBuffer.length; ++channel) {
		console.log(`Channel ${channel} has ${inputBuffer[channel].length} samples`)
		// Even though all channels should have the same ammount of samples...
	}
	return inputBuffer;
}

engine.addAudioCallback(processAudio)
```

### Alternatively, you can read/write samples to the sound card manually

```javascript
var engine = coreAudio.createNewAudioEngine()

// Grab the samples
var samples = engine.read()

// Silence the 0th channel
for(var sample = 0; sample < inputBuffer[0].length; ++ sample )
	samples[0][sample] = 0.0

// Send the sample array back to the sound card
engine.write(samples)
```

Important! Processing Thread
=====
When you are writing code inside of your audio callback, you are operating on
the processing thread of the application. This high priority environment means you
should try to think about performance as much as possible. Allocations and other
complex operations are possible, but dangerous.

IF YOU TAKE TOO LONG TO RETURN A BUFFER TO THE SOUND CARD, YOU WILL HAVE AUDIO DROPOUTS

The basic principle is that you should have everything ready to go before you enter
the processing function. Buffers, objects, and functions should be created in a constructor or static function outside of the audio callback whenever possible. The
examples in this readme are not necessarily good practice as far as performance is concerned.

The callback is only called if all buffers has been processed by the soundcard.

Audio Engine Options
=====

Option          | Type    | Default                  | Description
----------------|---------|--------------------------|----------------------------------------
sampleRate      | number  | 44100                    | The sample rate of the incoming audio.
sampleFormat    | string  | `sampleFormatFloat32`    | The sample format to use. Allowed `sampleFormatFloat32`, `sampleFormatInt32`, `sampleFormatInt24`, `sampleFormatInt16`, `sampleFormatInt8`, `sampleFormatUInt8`
framesPerBuffer | number  | 1024                     | The count of frames per buffer (per channel).
interleaved     | boolean | `false`                  | If set to `true` the samples are given as a two dimensional array `[channel][sample]`, otherwise the samples are given as a one dimensional array with alternating channel values.
inputChannels   | number  | 1                        | Count of input channels.
ouputChannels   | number  | 2                        | Count of ouput channels.
inputDevice     | number  | Default device id.        | Device id of the system's default input device.
outputDevice    | number  | Default output device id. | Device id of the system's default output device.

API
=====
First things first

```javascript
const coreAudio = require("node-core-audio")
```

Create and audio processing function

```javascript
function processAudio(samples) {
    // Just print the value of the first sample on the left channel
    console.log(samples[0][0])
    return samples
}
```

Initialize the audio engine and setup the processing loop

```javascript
var engine = coreAudio.createNewAudioEngine()
engine.addAudioCallback(processAudio)
```

General functionality

```javascript
/**
 * Returns whether the audio engine is active.
 * @return {boolean} Whether the audio engine is active.
 */
var active = engine.isActive()

/**
 * Updates the parameters and restarts the engine.
 * All keys from `getOptions()` are available.
 * @param {Object} options The audio engine options to set.
 */
engine.setOptions({
	inputChannels: 2
})

/**
 * Returns the current options for the engine.
 * @return {Object} The current audio engine options.
 */
var opts = engine.getOptions()

/**
 * Returns the current input buffer of the soundcard as array.
 * Note: This is a blocking call, don't take too long!
 * @return {Array} The current input buffer of the soundcard as array.
 */
var samples = engine.read()

/**
 * Writes an array of samples onto the soundcard.
 * Note: This is a blocking call, don't take too long!
 * @param {Array} samples An array of samples.
 * @return {boolean} If underflow occurred.
 */
var underflow = engine.write([0.42, ...])

/**
 * Returns the name of a given device id.
 * @param {number} device_id The device id to get the name of.
 * @return {string} The device name for the given id.
 */
var deviceName = engine.getDeviceName(0)

/**
 * Returns the total number of audio devices.
 * @return {number} The count of audio devices.
 */
var totalDevices = engine.getNumDevices()
```

Known Issues / TODO
=====

* Add FFTW to C++ extension, so you can get fast FFT's from javascript, and also register for the FFT of incoming audio, rather than the audio itself
* Add support for streaming audio over sockets


License
=====
MIT - See LICENSE file.

Copyright Mike Vegeto, 2013
