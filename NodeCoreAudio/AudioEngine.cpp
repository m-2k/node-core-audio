//////////////////////////////////////////////////////////////////////////////
//
// AudioEngine.cpp: Core audio functionality
//
//////////////////////////////////////////////////////////////////////////////
// DevStudio!namespace Audio DevStudio!class AudioEngine


#include "AudioEngine.h"
#include "portaudio.h"
#include <v8.h>
#include <uv.h>
#include <node.h>
#include <node_object_wrap.h>
#include <cstdlib>

#if !(defined( __WINDOWS__ ) || defined( _WIN32 ))
	#include <unistd.h>
	#include <cstring>
#endif

using namespace v8; using namespace std;

Persistent<Function> Audio::AudioEngine::constructor;

static void do_work( void* arg ) {
	Audio::AudioEngine* pEngine = (Audio::AudioEngine*)arg;
	pEngine->RunAudioLoop();
}

//////////////////////////////////////////////////////////////////////////////
/*! Initialize */
Audio::AudioEngine::AudioEngine( Local<Object> options ) :
	m_uSampleSize(4),
	m_pIsolate(Isolate::GetCurrent()),
	m_pLocker(new Locker(Isolate::GetCurrent())) {

	PaError openStreamErr;
	m_pPaStream = NULL;

	//set defaults
	m_uInputChannels = 2;
	m_uOutputChannels = 2;
	m_uSamplesPerBuffer = DEFAULT_FRAMES_PER_BUFFER;
	m_uNumCachedOutputSamples = 0;
	m_uSampleFormat = 1;
	m_uNumBuffers = DEFAULT_NUM_BUFFERS;
	m_bInterleaved = false;
	m_bReadMicrophone = true;

	m_uSampleRate = DEFAULT_SAMPLE_RATE;

	m_cachedInputSampleBlock = NULL;
	m_cachedOutputSampleBlock = NULL;
	m_cachedOutputSampleBlockForWriting = NULL;

	m_bInputOverflowed = 0;
	m_bOutputUnderflowed = 0;

	m_uCurrentWriteBuffer = 0;
	m_uCurrentReadBuffer = 0;

	// Create V8 objects to hold our buffers
	m_hInputBuffer = Nan::New<Array>( m_uInputChannels );
	for( int iChannel=0; iChannel<m_uInputChannels; iChannel++ )
		m_hInputBuffer->Set( iChannel, Nan::New<Array>(m_uSamplesPerBuffer) );

	// Initialize our audio core
	PaError initErr = Pa_Initialize();

	m_uInputDevice = Pa_GetDefaultInputDevice();
	if( m_uInputDevice == paNoDevice ) {
		Nan::ThrowTypeError("Error: No default input device");
	};

	m_uOutputDevice = Pa_GetDefaultOutputDevice();
	if( m_uOutputDevice == paNoDevice ) {
		Nan::ThrowTypeError("Error: No default output device");
	}

	if( initErr != paNoError )
		Nan::ThrowTypeError("Failed to initialize audio engine");

	// Set fft defaults
	fft_callback.Reset();
	fft_window_size = 1024;
	fft_overlap_size = 0.0;
	fft_window_function_type = Blackman;

	fft_resources_allocated = false;

	applyOptions( options );	

	// Open an audio I/O stream.
	openStreamErr = Pa_OpenStream(  &m_pPaStream,
									&m_inputParams,
									&m_outputParams,
									m_uSampleRate,
									m_uSamplesPerBuffer,
									paClipOff,
									NULL,
									NULL );

	if( openStreamErr != paNoError )
		Nan::ThrowTypeError("Failed to open audio stream");

	// Start the audio stream
	PaError startStreamErr = Pa_StartStream( m_pPaStream );

	if( startStreamErr != paNoError )
		Nan::ThrowTypeError("Failed to start audio stream");

	uv_mutex_init( &m_mutex );
	uv_thread_create( &ptStreamThread, do_work, (void*)this );

} // end Constructor

/**
 * Allocates the fft resources for input stream fft computations.
 */
void Audio::AudioEngine::setupFFT() {
	if (fft_resources_allocated) freeFFT();

	//printf("Allocating fft stuff\n");
	fft_window_function = new WindowFunction(fft_window_function_type, fft_window_size);

	fft_windows = new double*[m_uInputChannels];
	fft_window_idx = new int[m_uInputChannels];

	fft_results = (fftw_complex**) malloc(m_uInputChannels * sizeof(fftw_complex*));
	fft_plans = (fftw_plan*) malloc(m_uInputChannels * sizeof(fftw_plan));
	
	for (int iChannel = 0; iChannel < m_uInputChannels; ++iChannel) {
		fft_windows[iChannel] = new double[fft_window_size];
		fft_window_idx[iChannel] = 0;
		fft_results[iChannel] = (fftw_complex*) fftw_malloc(((fft_window_size / 2) + 1) * sizeof(fftw_complex));
		fft_plans[iChannel] = fftw_plan_dft_r2c_1d(fft_window_size, fft_windows[iChannel], fft_results[iChannel], FFTW_ESTIMATE);
	}
	fft_resources_allocated = true;
}


//////////////////////////////////////////////////////////////////////////////
/*! Gets options */
void Audio::AudioEngine::getOptions(const Nan::FunctionCallbackInfo<v8::Value>& info) {
	Local<Object> options = Nan::New<Object>();

	AudioEngine* pEngine = AudioEngine::Unwrap<AudioEngine>( info.This() );

	Nan::Set( options, Nan::New<String>("inputChannels").ToLocalChecked(), Nan::New<Number>(pEngine->m_uInputChannels) );
	Nan::Set( options, Nan::New<String>("outputChannels").ToLocalChecked(), Nan::New<Number>(pEngine->m_uOutputChannels) );

	Nan::Set( options, Nan::New<String>("inputDevice").ToLocalChecked(), Nan::New<Number>(pEngine->m_uInputDevice) );
	Nan::Set( options, Nan::New<String>("outputDevice").ToLocalChecked(), Nan::New<Number>(pEngine->m_uOutputDevice) );

	Nan::Set( options, Nan::New<String>("sampleRate").ToLocalChecked(), Nan::New<Number>(pEngine->m_uSampleRate) );
	Nan::Set( options, Nan::New<String>("sampleFormat").ToLocalChecked(), Nan::New<Number>(pEngine->m_uSampleFormat) );
	Nan::Set( options, Nan::New<String>("framesPerBuffer").ToLocalChecked(), Nan::New<Number>(pEngine->m_uSamplesPerBuffer) );
	Nan::Set( options, Nan::New<String>("numBuffers").ToLocalChecked(), Nan::New<Number>(pEngine->m_uNumBuffers) );
	Nan::Set( options, Nan::New<String>("interleaved").ToLocalChecked(), Nan::New<Boolean>(pEngine->m_bInterleaved) );
	Nan::Set( options, Nan::New<String>("useMicrophone").ToLocalChecked(), Nan::New<Boolean>(pEngine->m_bReadMicrophone) );

	info.GetReturnValue().Set(options);
} // end GetOptions


//////////////////////////////////////////////////////////////////////////////
/*! Set options, restarts audio stream */
void Audio::AudioEngine::setOptions(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;
	Local<Object> options;

	if( info.Length() > 0 ) {
		if( !info[0]->IsObject() ) {
			return Nan::ThrowTypeError("First argument should be an object.");
		}

		options = Local<Object>::Cast( info[0] );
	} else {
        return Nan::ThrowTypeError("First argument does not exist.");
	}

	AudioEngine* pEngine = AudioEngine::Unwrap<AudioEngine>( info.This() );
	pEngine->freeFFT();
	pEngine->applyOptions( options );
	pEngine->setupFFT();

    	info.GetReturnValue().SetUndefined();
} // end AudioEngine::SetOptions()x


//////////////////////////////////////////////////////////////////////////////
/*! Sets the given options and restarts the audio stream if necessary */
void Audio::AudioEngine::applyOptions( Local<Object> options ) {
	unsigned int oldBufferCount = m_uNumBuffers;
	if(Nan::HasOwnProperty(options, Nan::New<String>("inputDevice").ToLocalChecked()).FromMaybe(false) )
		m_uInputDevice = Nan::To<int>(Nan::Get(options, Nan::New<String>("inputDevice").ToLocalChecked()).ToLocalChecked()).FromJust();
	if(Nan::HasOwnProperty(options, Nan::New<String>("outputDevice").ToLocalChecked()).FromMaybe(false) )
		m_uOutputDevice = Nan::To<int>(Nan::Get(options, Nan::New<String>("outputDevice").ToLocalChecked()).ToLocalChecked()).FromJust();
	if(Nan::HasOwnProperty(options, Nan::New<String>("inputChannels").ToLocalChecked()).FromMaybe(false) )
		m_uInputChannels = Nan::To<int>(Nan::Get(options, Nan::New<String>("inputChannels").ToLocalChecked()).ToLocalChecked()).FromJust();
	if(Nan::HasOwnProperty(options, Nan::New<String>("outputChannels").ToLocalChecked()).FromMaybe(false) )
	  m_uOutputChannels = Nan::To<int>(Nan::Get(options, Nan::New<String>("outputChannels").ToLocalChecked()).ToLocalChecked()).FromJust();
	if(Nan::HasOwnProperty(options, Nan::New<String>("framesPerBuffer").ToLocalChecked()).FromMaybe(false) )
		m_uSamplesPerBuffer = Nan::To<int>(Nan::Get(options, Nan::New<String>("framesPerBuffer").ToLocalChecked()).ToLocalChecked()).FromJust();
	if (Nan::HasOwnProperty(options, Nan::New<String>("numBuffers").ToLocalChecked()).FromMaybe(false) )
		m_uNumBuffers = Nan::To<int>(Nan::Get(options, Nan::New<String>("numBuffers").ToLocalChecked()).ToLocalChecked()).FromJust();
	if(Nan::HasOwnProperty(options, Nan::New<String>("interleaved").ToLocalChecked()).FromMaybe(false) )
		m_bInterleaved = Nan::To<bool>(Nan::Get(options, Nan::New<String>("interleaved").ToLocalChecked()).ToLocalChecked()).FromJust();
	if (Nan::HasOwnProperty(options, Nan::New<String>("useMicrophone").ToLocalChecked()).FromMaybe(false) )
		m_bReadMicrophone = Nan::To<bool>(Nan::Get(options, Nan::New<String>("useMicrophone").ToLocalChecked()).ToLocalChecked()).FromJust();
	if (Nan::HasOwnProperty(options, Nan::New<String>("sampleRate").ToLocalChecked()).FromMaybe(false) )
		m_uSampleRate = Nan::To<int>(Nan::Get(options, Nan::New<String>("sampleRate").ToLocalChecked()).ToLocalChecked()).FromJust();
	if(Nan::HasOwnProperty(options, Nan::New<String>("sampleFormat").ToLocalChecked()).FromMaybe(false) ) {
		switch(Nan::To<int>(Nan::Get(options, Nan::New<String>("sampleFormat").ToLocalChecked()).ToLocalChecked()).FromJust()){
			case 0x01: m_uSampleFormat = paFloat32; m_uSampleSize = 4; break;
			case 0x02: m_uSampleFormat = paInt32; m_uSampleSize = 4; break;
			case 0x04: m_uSampleFormat = paInt24; m_uSampleSize = 3; break;
			case 0x08: m_uSampleFormat = paInt16; m_uSampleSize = 2; break;
			case 0x10: m_uSampleFormat = paInt8; m_uSampleSize = 1; break;
			case 0x20: m_uSampleFormat = paUInt8; m_uSampleSize = 1; break;
		}
	}

	
	// Free the resources that where allocated before
	freeFFT();

	// FFT options
	if (Nan::HasOwnProperty(options, Nan::New<String>("fftWindowSize").ToLocalChecked()).FromMaybe(false) ) {
		int desiredWindowSize = Nan::To<int>(Nan::Get(options, Nan::New<String>("fftWindowSize").ToLocalChecked()).ToLocalChecked()).FromJust();
		if (desiredWindowSize < 256 || desiredWindowSize > 4096 || (desiredWindowSize & (desiredWindowSize - 1)) != 0) {
			printf("Invalid fft window size %i (fft window size must be a power of two between 256 and 4096).\n", desiredWindowSize);
		} else {
			fft_window_size = desiredWindowSize;
		}
	}
	if (Nan::HasOwnProperty(options, Nan::New<String>("fftOverlapSize").ToLocalChecked()).FromMaybe(false) ) {
		float desiredOverlapSize = (float)Nan::To<double>(Nan::Get(options, Nan::New<String>("fftOverlapSize").ToLocalChecked()).ToLocalChecked()).FromJust();
		if (desiredOverlapSize >= 1 || desiredOverlapSize < 0) {
			printf("Invalid fft overlap size %f (fft overlap size must be a number >= 0 and < 1)\n", desiredOverlapSize);
		} else {
			fft_overlap_size = desiredOverlapSize;
		}
	}
	if (Nan::HasOwnProperty(options, Nan::New<String>("fftWindowFunction").ToLocalChecked()).FromMaybe(false) ) {
		Local<String> _str = Nan::To<String>(Nan::Get(options, Nan::New<String>("fftWindowFunction").ToLocalChecked()).ToLocalChecked()).ToLocalChecked();
		char* desiredWindowFunction = (char*)(*String::Utf8Value(_str));
		
		
			if (strcmp(desiredWindowFunction, "Square") == 0)				fft_window_function_type = Square;
			else if (strcmp(desiredWindowFunction, "VonHann") == 0)			fft_window_function_type = VonHann;
			else if (strcmp(desiredWindowFunction, "Hamming") == 0)			fft_window_function_type = Hamming;
			else if (strcmp(desiredWindowFunction, "Blackman") == 0) 		fft_window_function_type = Blackman;
			else if (strcmp(desiredWindowFunction, "BlackmanHarris") == 0)	fft_window_function_type = BlackmanHarris;
			else if (strcmp(desiredWindowFunction, "BlackmanNuttall") == 0)	fft_window_function_type = BlackmanNuttall;
			else if (strcmp(desiredWindowFunction, "FlatTop") == 0)			fft_window_function_type = FlatTop;
			else printf("Unknown window function %s.\n", desiredWindowFunction);
	}


	// Setup our input and output parameters
	m_inputParams.device = m_uInputDevice;
	m_inputParams.channelCount = m_uInputChannels;
	m_inputParams.sampleFormat = m_uSampleFormat;
	m_inputParams.suggestedLatency = Pa_GetDeviceInfo(m_inputParams.device)->defaultLowInputLatency;
	m_inputParams.hostApiSpecificStreamInfo = NULL;

	m_outputParams.device = m_uOutputDevice;
	m_outputParams.channelCount = m_uOutputChannels;
	m_outputParams.sampleFormat = m_uSampleFormat;
	m_outputParams.suggestedLatency = Pa_GetDeviceInfo(m_outputParams.device)->defaultLowOutputLatency;
	m_outputParams.hostApiSpecificStreamInfo = NULL;

	//printf("Interleaved %d\n", m_bInterleaved);

	// Clear out our temp buffer blocks
	if( m_cachedInputSampleBlock != NULL )
		free( m_cachedInputSampleBlock );
	if( m_cachedOutputSampleBlock != NULL ) {
		for (unsigned int i = 0; i < oldBufferCount; i++) {
			if ( m_cachedOutputSampleBlock[i] != NULL ) {
				free(m_cachedOutputSampleBlock[i]);
			}
		}
		free( m_cachedOutputSampleBlock );
	}
	if ( m_cachedOutputSampleBlockForWriting != NULL)
		free (m_cachedOutputSampleBlockForWriting);

	// Allocate some new space for our temp buffer blocks
	m_cachedInputSampleBlock = (char*)malloc( m_uSamplesPerBuffer * m_uInputChannels * m_uSampleSize );
	m_cachedOutputSampleBlock = (char**)malloc( sizeof(char*) * m_uNumBuffers);
	for (int i = 0; i < m_uNumBuffers; i++) {
		m_cachedOutputSampleBlock[i] = (char*)malloc( m_uSamplesPerBuffer * m_uOutputChannels * m_uSampleSize );
	}
	m_cachedOutputSampleBlockForWriting = (char*)malloc( m_uSamplesPerBuffer * m_uOutputChannels * m_uSampleSize );
	m_uNumCachedOutputSamples = (unsigned int*)calloc( sizeof(unsigned int), m_uNumBuffers );

	// Reallocate the fft resources for the new settings
	setupFFT();

	if( m_pPaStream != NULL && Pa_IsStreamActive(m_pPaStream) )
		restartStream();

} // end AudioEngine::applyOptions()

//////////////////////////////////////////////////////////////////////////////
/*! Returns a v8 array filled with input samples */
Local<Array> Audio::AudioEngine::getInputBuffer() {
  Nan::EscapableHandleScope scope;

	if( m_bInterleaved ) {
		m_hInputBuffer = Nan::New<Array>( m_uInputChannels * m_uSamplesPerBuffer );

		for( int iSample=0; iSample<m_uSamplesPerBuffer*m_uInputChannels; ++iSample ) {
			m_hInputBuffer->Set( iSample, getSample(iSample) );
		}
	} else {
		// The NonInterleaved flag will never be set, therefore 
		// we create a non interleaved array [channel][sample]
		m_hInputBuffer = Nan::New<Array>(m_uInputChannels);
		for(int iChannel = 0; iChannel < m_uInputChannels; ++iChannel) {
			auto tempBuffer = Local<Array>(Nan::New<Array>(m_uSamplesPerBuffer));
			int rSample = 0;
			for (int iSample = iChannel; iSample < m_uSamplesPerBuffer; iSample += m_uInputChannels) {
				tempBuffer->Set(rSample, getSample(iSample));
				++rSample;
			}
			m_hInputBuffer->Set(iChannel, tempBuffer);
		}

		//m_hInputBuffer = Nan::New<Array>( m_uInputChannels );
		/*for( int iChannel=0; iChannel<m_uInputChannels; iChannel++ ) {
			auto tempBuffer = Local<Array>( Nan::New<Array>(m_uSamplesPerBuffer) );

			for( int iSample=0; iSample<m_uSamplesPerBuffer; iSample++ ) {
				tempBuffer->Set( iSample, getSample(iSample) );
			}

			m_hInputBuffer->Set( iChannel, tempBuffer );
		}*/
	}

	return scope.Escape( m_hInputBuffer );
} // end AudioEngine::getInputBuffer()


//////////////////////////////////////////////////////////////////////////////
/*! Returns a sound card sample converted to a v8 Number */
Handle<Number> Audio::AudioEngine::getSample( int position ) {
  Nan::EscapableHandleScope scope;

	Local<Number> sample;

	switch( m_uSampleFormat ) {
	case paFloat32: {
		float fValue = ((float*)m_cachedInputSampleBlock)[position];

		sample = Nan::New<Number>( fValue );
		break;
					}
	case paInt32:
		sample = Nan::New<Integer>( ((int*)m_cachedInputSampleBlock)[position] );
		break;

	case paInt24:
		sample = Nan::New<Integer>(
				(m_cachedInputSampleBlock[3*position + 0] << 16)
			  + (m_cachedInputSampleBlock[3*position + 1] << 8)
			  + (m_cachedInputSampleBlock[3*position + 2]) );
		break;

	case paInt16:
		sample = Nan::New<Integer>( ((int16_t*)m_cachedInputSampleBlock)[position] );
		break;

	default:
		sample = Nan::New<Integer>( m_cachedInputSampleBlock[position]*-1 );
		break;
	}

    return scope.Escape(sample);
} // end AudioEngine::getSample()


/**
 * Sets the fft callback function.
 */
void Audio::AudioEngine::setFFTCallback(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;

	AudioEngine* pEngine = AudioEngine::Unwrap<AudioEngine>( info.This() );

	if (info.Length() > 1 || !info[0]->IsFunction()){
        return Nan::ThrowTypeError("First argument must be a callback function.");
	}

	//pEngine->fft_callback(Local<Function>::Cast(info[0]));
	//pEngine->fft_callback.SetFunction(Local<Function>::Cast(info[0]));
	pEngine->fft_callback.Reset(Local<Function>::Cast(info[0]));
}

/**
 * Returns a sound card sample converted to a double.
 */
double Audio::AudioEngine::getDoubleSample( int position ) {
	switch( m_uSampleFormat ) {
		case paFloat32:
			return (double)((float*)m_cachedInputSampleBlock)[position];
		case paInt32:
			return (double)((int*)m_cachedInputSampleBlock)[position];
		case paInt24:
			return (double)((m_cachedInputSampleBlock[3*position + 0] << 16)
							+ (m_cachedInputSampleBlock[3*position + 1] << 8)
							+ (m_cachedInputSampleBlock[3*position + 2]));
		case paInt16:
			return (double)((int16_t*)m_cachedInputSampleBlock)[position];
	}
	return (double)m_cachedInputSampleBlock[position]*-1;
} // end AudioEngine::getSample()

/**
 * Fills the fft windows and executes them when they are full.
 */
void Audio::AudioEngine::fillWindows() {
	if (fft_resources_allocated == false || fft_callback.IsEmpty()) return;

	// Fill the channel fft windows
	for (int iSample = 0; iSample < m_uSamplesPerBuffer * m_uInputChannels; ++iSample) {
		int iChannel = iSample % m_uInputChannels;
		double sample = getDoubleSample(iSample + iChannel);
		fft_windows[iChannel][fft_window_idx[iChannel]] = sample * fft_window_function->at(fft_window_idx[iChannel]);
		
		// Check if the fft window for this channel is full and the fft plan can be executed
		++fft_window_idx[iChannel];
		if (fft_window_idx[iChannel] >= fft_window_size) {
			runPlan(iChannel);
			fft_window_idx[iChannel] = 0;
			
			// Add the overlapping last values to the window again
			if (fft_overlap_size > 0.0 && fft_overlap_size <= 1.0) {
				int overlapping_samples = (int) floor(fft_overlap_size * (float)fft_window_size);
				for (int iOverlap = fft_window_size - overlapping_samples; iOverlap < fft_window_size; ++iOverlap) {
					fft_windows[iChannel][fft_window_idx[iChannel]] = fft_windows[iChannel][iOverlap];
					++fft_window_idx[iChannel];
				}
			}

			// Create an array of the result and normalize it
			int result_length = (fft_window_size / 2) + 1;

			Local<Array> results = Nan::New<Array>(result_length);
			for (int i = 0; i < result_length; ++i) {
				//results[i] = fft_results[iChannel][i][0];
				double db = 10.0f * log(fft_results[iChannel][i][0]*fft_results[iChannel][i][0] + fft_results[iChannel][i][1]*fft_results[iChannel][i][1]) / log(10.0);
				results->Set(i, Nan::New<Number>(db));
			}
			
			v8::Local<v8::Value> argv[2] = {Nan::New<Number>(iChannel), results};
			fft_callback.Call(2, argv);
		}
	}
}

/**
 * Executes the plan for a specific channel.
 */
void Audio::AudioEngine::runPlan(int iChannel) {
	fftw_execute(fft_plans[iChannel]);
}


//////////////////////////////////////////////////////////////////////////////
/*! Sets a sample in the queued output buffer */
void Audio::AudioEngine::setSample( int position, Handle<Value> sample ) {
	int temp;
	switch( m_uSampleFormat ) {
	case paFloat32:
		((float *)m_cachedOutputSampleBlock[m_uCurrentWriteBuffer])[position] = (float)sample->NumberValue();
		break;

	case paInt32:
		((int *)m_cachedOutputSampleBlock[m_uCurrentWriteBuffer])[position] = (int)sample->NumberValue();
		break;

	case paInt24:
		temp = (int)sample->NumberValue();
		m_cachedOutputSampleBlock[m_uCurrentWriteBuffer][3*position + 0] = temp >> 16;
		m_cachedOutputSampleBlock[m_uCurrentWriteBuffer][3*position + 1] = temp >> 8;
		m_cachedOutputSampleBlock[m_uCurrentWriteBuffer][3*position + 2] = temp;
		break;

	case paInt16:
		((int16_t *)m_cachedOutputSampleBlock[m_uCurrentWriteBuffer])[position] = (int16_t)sample->NumberValue();
		break;

	default:
		m_cachedOutputSampleBlock[m_uCurrentWriteBuffer][position] = (char)sample->NumberValue();
		break;
	}
} // end AudioEngine::setSample()


//////////////////////////////////////////////////////////////////////////////
/*! Queues up an array to be sent to the sound card */
void Audio::AudioEngine::queueOutputBuffer( Handle<Array> result ) {
	// Reset our record of the number of cached output samples
	m_uNumCachedOutputSamples[m_uCurrentWriteBuffer] = 0;

	if( m_bInterleaved ) {
		for( int iSample=0; iSample<m_uSamplesPerBuffer*m_uOutputChannels; ++iSample )
			setSample( iSample, result->Get(iSample) );

		m_uNumCachedOutputSamples[m_uCurrentWriteBuffer] = result->Length()/m_uOutputChannels;

	} else {
		Local<Array> item;
		for(int iChannel = 0; iChannel < m_uOutputChannels; ++iChannel) {
			item = Local<Array>::Cast(result->Get(iChannel));
			if (!item->IsArray()) {
				char errStr [50];
				sprintf(errStr, "Output buffer not properly setup, %ith channel is not an array", iChannel);
				Nan::ThrowTypeError(errStr);
				return;
			}
			if (item->Length() > m_uNumCachedOutputSamples[m_uCurrentWriteBuffer])
				m_uNumCachedOutputSamples[m_uCurrentWriteBuffer] = item->Length();
			int rSample = iChannel;
			for(int iSample = 0; iSample < m_uSamplesPerBuffer; ++iSample) {
				rSample += m_uOutputChannels;
				setSample(rSample, item->Get(iSample));
			}
		}
		/*for( int iChannel=0; iChannel<m_uOutputChannels; ++iChannel ) {
			for( int iSample=0; iSample<m_uSamplesPerBuffer; ++iSample ) {

				item = Local<Array>::Cast( result->Get(iChannel) );
				if( item->IsArray() ) {
					if( item->Length() > m_uNumCachedOutputSamples[m_uCurrentWriteBuffer] )
						m_uNumCachedOutputSamples[m_uCurrentWriteBuffer] = item->Length();

					setSample( iSample, item->Get(iSample) );
				}
			} // end for each sample
		} // end for each channel*/
	}
	m_uCurrentWriteBuffer = (m_uCurrentWriteBuffer + 1)%m_uNumBuffers;
} // end AudioEngine::queueOutputBuffer()


//////////////////////////////////////////////////////////////////////////////
/*! Run the main blocking audio loop */
void Audio::AudioEngine::RunAudioLoop(){

	PaError error = 0;

	assert( m_pPaStream );

	while( true ) {
		m_bInputOverflowed = (error != paNoError);

		int numSamples = 0;
		uv_mutex_lock( &m_mutex );
		if (m_uNumCachedOutputSamples[m_uCurrentReadBuffer] > 0) {
			memcpy(m_cachedOutputSampleBlockForWriting, m_cachedOutputSampleBlock[m_uCurrentReadBuffer], m_uNumCachedOutputSamples[m_uCurrentReadBuffer] * m_uOutputChannels * m_uSampleSize);
			numSamples = m_uNumCachedOutputSamples[m_uCurrentReadBuffer];
			m_uNumCachedOutputSamples[m_uCurrentReadBuffer] = 0;
			m_uCurrentReadBuffer = (m_uCurrentReadBuffer + 1)%m_uNumBuffers;
		}
		uv_mutex_unlock( &m_mutex );

		if( numSamples > 0 ) {
			error = Pa_WriteStream( m_pPaStream, m_cachedOutputSampleBlockForWriting, numSamples );
			m_bOutputUnderflowed = (error != paNoError);
		} else {
			m_bOutputUnderflowed = true;
#if defined( __WINDOWS__ ) || defined( _WIN32 )
			Sleep(1);
#else
			sleep(1);
#endif
		}
	}
} // end AudioEngine::RunAudioLoop()


//////////////////////////////////////////////////////////////////////////////
/*! Initialize our node object */
NAN_MODULE_INIT(Audio::AudioEngine::Init) {

	// Prepare constructor template
	Local<FunctionTemplate> functionTemplate = Nan::New<FunctionTemplate> (Audio::AudioEngine::New );
	functionTemplate->SetClassName( Nan::New<String>("AudioEngine").ToLocalChecked() );
	functionTemplate->InstanceTemplate()->SetInternalFieldCount( 1 );


    //Local<FunctionTemplate> constructorHandle = Nan::New(constructor);
    //target->Set(Nan::New<String>("AudioEngine"), functionTemplate->GetFunction());

    // Get
	//functionTemplate->PrototypeTemplate()->Set( Nan::New<String>("isActive"), Nan::New<FunctionTemplate>(Audio::AudioEngine::isActive)->GetFunction() );
	//functionTemplate->PrototypeTemplate()->Set( Nan::New<String>("getDeviceName"), Nan::New<FunctionTemplate>(Audio::AudioEngine::getDeviceName)->GetFunction() );
	//functionTemplate->PrototypeTemplate()->Set( Nan::New<String>("getNumDevices"), Nan::New<FunctionTemplate>(Audio::AudioEngine::getNumDevices)->GetFunction() );
    Nan::SetPrototypeMethod(functionTemplate, "isActive", Audio::AudioEngine::isActive);
    Nan::SetPrototypeMethod(functionTemplate, "getDeviceName", Audio::AudioEngine::getDeviceName);
    Nan::SetPrototypeMethod(functionTemplate, "getDeviceMaxInputChannels", Audio::AudioEngine::getDeviceMaxInputChannels);
    Nan::SetPrototypeMethod(functionTemplate, "getDeviceMaxOutputChannels", Audio::AudioEngine::getDeviceMaxOutputChannels);
    Nan::SetPrototypeMethod(functionTemplate, "getNumDevices", Audio::AudioEngine::getNumDevices);

	// Set
	//functionTemplate->PrototypeTemplate()->Set( Nan::New<String>("setOptions"), Nan::New<FunctionTemplate>(Audio::AudioEngine::setOptions)->GetFunction() );
	//functionTemplate->PrototypeTemplate()->Set( Nan::New<String>("getOptions"), Nan::New<FunctionTemplate>(Audio::AudioEngine::getOptions)->GetFunction() );
	//functionTemplate->PrototypeTemplate()->Set( Nan::New<String>("write"), Nan::New<FunctionTemplate>(Audio::AudioEngine::write)->GetFunction() );
	//functionTemplate->PrototypeTemplate()->Set( Nan::New<String>("read"), Nan::New<FunctionTemplate>(Audio::AudioEngine::read)->GetFunction() );
	//functionTemplate->PrototypeTemplate()->Set( Nan::New<String>("isBufferEmpty"), Nan::New<FunctionTemplate>(Audio::AudioEngine::isBufferEmpty)->GetFunction() );
    Nan::SetPrototypeMethod(functionTemplate, "setOptions", Audio::AudioEngine::setOptions);
    Nan::SetPrototypeMethod(functionTemplate, "getOptions", Audio::AudioEngine::getOptions);
    Nan::SetPrototypeMethod(functionTemplate, "write", Audio::AudioEngine::write);
    Nan::SetPrototypeMethod(functionTemplate, "read", Audio::AudioEngine::read);
    Nan::SetPrototypeMethod(functionTemplate, "isBufferEmpty", Audio::AudioEngine::isBufferEmpty);

    Nan::SetPrototypeMethod(functionTemplate, "setFFTCallback", Audio::AudioEngine::setFFTCallback);

	//constructor = Persistent<Function>::New( functionTemplate->GetFunction() );
    //Local<FunctionTemplate> tpl = Nan::New<FunctionTemplate>(EOLFinder::New);
//    NanAssignPersistent(constructor, functionTemplate->GetFunction());
    constructor.Reset(Isolate::GetCurrent(), functionTemplate->GetFunction());
} // end AudioEngine::Init()


//////////////////////////////////////////////////////////////////////////////
/*! Create a new instance of the audio engine */
//v8::Handle<v8::Value> Audio::AudioEngine::NewInstance(const v8::Arguments& info) {
void Audio::AudioEngine::NewInstance(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	unsigned argc = info.Length();

	if( argc > 2 )
		argc = 2;

	Handle<Value>* argv = new Handle<Value>[argc];

	argv[0] = info[0];
	if( argc > 1 )
		argv[1] = info[1];

	//Local<Object> instance = constructor->NewInstance( argc, argv );
	Local<Object> instance = Nan::New(constructor)->NewInstance(argc, argv);
	//Local<Object> instance = constructor->NewInstance(argc, argv);

	info.GetReturnValue().Set( instance );
} // end AudioEngine::NewInstance()

//////////////////////////////////////////////////////////////////////////////
/*! Create a v8 object */
void Audio::AudioEngine::New(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	Local<Object> options;

	if( info.Length() > 0 ) {
		if( !info[0]->IsObject() )
            return Nan::ThrowTypeError("First argument must be an object.");
		else
			options = Local<Object>::Cast( info[0] );
	} else {
		options = Nan::New<Object>();
	}

	AudioEngine* pEngine = new AudioEngine( options );
	pEngine->Wrap( info.This() );

	info.GetReturnValue().Set( info.This() );
} // end AudioEngine::New()


//////////////////////////////////////////////////////////////////////////////
/*! Write samples to the current audio device */
void Audio::AudioEngine::write(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	AudioEngine* pEngine = AudioEngine::Unwrap<AudioEngine>( info.This() );

	if (info.Length() > 1 || !info[0]->IsArray()){
        return Nan::ThrowTypeError("First argument should be an array.");
	}

	uv_mutex_lock( &pEngine->m_mutex );
	pEngine->queueOutputBuffer( Local<Array>::Cast(info[0]) );
	uv_mutex_unlock( &pEngine->m_mutex );

	info.GetReturnValue().Set( false );
} // end AudioEngine::Write()

//////////////////////////////////////////////////////////////////////////////
/*! Checks if the current audio buffer has been fed to Port Audio */
void Audio::AudioEngine::isBufferEmpty(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	AudioEngine* pEngine = AudioEngine::Unwrap<AudioEngine>( info.This() );

	uv_mutex_lock( &pEngine->m_mutex );
	Local<Boolean> isEmpty = Nan::New<Boolean>(pEngine->m_uNumCachedOutputSamples[pEngine->m_uCurrentWriteBuffer] == 0);
	uv_mutex_unlock( &pEngine->m_mutex );
	info.GetReturnValue().Set( isEmpty );
} // end AudioEngine::isBufferEmpty()

//////////////////////////////////////////////////////////////////////////////
/*! Read samples from the current audio device */
void Audio::AudioEngine::read(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	AudioEngine* pEngine = AudioEngine::Unwrap<AudioEngine>( info.This() );

	if (pEngine->m_bReadMicrophone) {
		Pa_ReadStream( pEngine->m_pPaStream, pEngine->m_cachedInputSampleBlock, pEngine->m_uSamplesPerBuffer );
	}

	Local<Array> input = pEngine->getInputBuffer();

	// Add the samples to the fft windows and execute them if possible
	pEngine->fillWindows();

	info.GetReturnValue().Set( input );
} // end AudioEngine::Read()


//////////////////////////////////////////////////////////////////////////////
/*! Returns whether the PortAudio stream is active */
void Audio::AudioEngine::isActive(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	AudioEngine* pEngine = AudioEngine::Unwrap<AudioEngine>( info.This() );

	if( Pa_IsStreamActive(pEngine->m_pPaStream) )
		info.GetReturnValue().Set( Nan::New<Boolean>(true) );
	else
		info.GetReturnValue().Set( Nan::New<Boolean>(false) );
} // end AudioEngine::IsActive()


//////////////////////////////////////////////////////////////////////////////
/*! Get the name of an audio device with a given ID number */
void Audio::AudioEngine::getDeviceName(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	if( !info[0]->IsNumber() ) {
		return  Nan::ThrowTypeError("getDeviceName() requires a device index");
	}

	Local<Number> deviceIndex = Local<Number>::Cast( info[0] );

	const PaDeviceInfo* pDeviceInfo = Pa_GetDeviceInfo( (PaDeviceIndex)deviceIndex->NumberValue() );

	info.GetReturnValue().Set( Nan::New<String>(pDeviceInfo->name).ToLocalChecked() );
} // end AudioEngine::GetDeviceName()

//////////////////////////////////////////////////////////////////////////////
/*! Get the name of an audio device with a given ID number */
void Audio::AudioEngine::getDeviceMaxInputChannels(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	if( !info[0]->IsNumber() ) {
		return  Nan::ThrowTypeError("getDeviceMaxInputChannels() requires a device index");
	}

	Local<Number> deviceIndex = Local<Number>::Cast( info[0] );

	const PaDeviceInfo* pDeviceInfo = Pa_GetDeviceInfo( (PaDeviceIndex)deviceIndex->NumberValue() );

	info.GetReturnValue().Set( Nan::New<Number>(pDeviceInfo->maxInputChannels) );
} // end AudioEngine::getDeviceMaxInputChannels()

//////////////////////////////////////////////////////////////////////////////
/*! Get the name of an audio device with a given ID number */
void Audio::AudioEngine::getDeviceMaxOutputChannels(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	if( !info[0]->IsNumber() ) {
		return  Nan::ThrowTypeError("getDeviceMaxOutputChannels() requires a device index");
	}

	Local<Number> deviceIndex = Local<Number>::Cast( info[0] );

	const PaDeviceInfo* pDeviceInfo = Pa_GetDeviceInfo( (PaDeviceIndex)deviceIndex->NumberValue() );

	info.GetReturnValue().Set( Nan::New<Number>(pDeviceInfo->maxOutputChannels) );
} // end AudioEngine::getDeviceMaxOutputChannels()



//////////////////////////////////////////////////////////////////////////////
/*! Get the number of available devices */
void Audio::AudioEngine::getNumDevices(const Nan::FunctionCallbackInfo<v8::Value>& info) {
    Nan::HandleScope scope;
	//HandleScope scope;

	int deviceCount = Pa_GetDeviceCount();

	info.GetReturnValue().Set( Nan::New<Number>(deviceCount) );
} // end AudioEngine::GetNumDevices()


//////////////////////////////////////////////////////////////////////////////
/*! Closes and reopens the PortAudio stream */
void Audio::AudioEngine::restartStream() {
	PaError error;

	// Stop the audio stream
	error = Pa_StopStream( m_pPaStream );

	if( error != paNoError )
		Nan::ThrowTypeError("Failed to stop audio stream");

	// Close the audio stream
	error = Pa_CloseStream( m_pPaStream );

	if( error != paNoError )
		Nan::ThrowTypeError("Failed to close audio stream");

	// Open an audio stream.``
	error = Pa_OpenStream(  &m_pPaStream,
		&m_inputParams,
		&m_outputParams,
		m_uSampleRate,
		m_uSamplesPerBuffer,
		paClipOff,              // We won't output out of range samples so don't bother clipping them
		NULL,
		NULL );

	if( error != paNoError ) {
		Nan::ThrowTypeError("Failed to open audio stream :(");
	}

	// Start the audio stream
	error = Pa_StartStream( m_pPaStream );

	if( error != paNoError )
		Nan::ThrowTypeError("Failed to start audio stream");

} // end AudioEngine::restartStream()


//////////////////////////////////////////////////////////////////////////////
/*! Wraps a handle into an object */
void Audio::AudioEngine::wrapObject( v8::Handle<v8::Object> object ) {
	ObjectWrap::Wrap( object );
} // end AudioEngine::wrapObject()

/**
 * Frees the fft resources	
 */
void Audio::AudioEngine::freeFFT() {
	// Check if the resources have been allocated
	if (fft_resources_allocated == false) return;
	//printf("Deallocating fft stuff\n");
    for (int iChannel = 0; iChannel < m_uInputChannels; ++iChannel) {
		delete[] fft_windows[iChannel];
		fftw_free(fft_results[iChannel]);
		fftw_destroy_plan(fft_plans[iChannel]);
	}
	delete[] fft_windows;
	delete[] fft_window_idx;
	free(fft_results);
	free(fft_plans);
	delete fft_window_function;
	fft_resources_allocated = false;
}

//////////////////////////////////////////////////////////////////////////////
/*! OOL Destructor */
Audio::AudioEngine::~AudioEngine() {
	// Kill PortAudio
	PaError err = Pa_Terminate();

	if( err != paNoError )
		fprintf( stderr, "PortAudio error: %s\n", Pa_GetErrorText( err ) );

	// Free fft resources
    freeFFT();
} // end AudioEngine::~AudioEngine()
