// Create a new instance of node-core-audio
const coreAudio = require('../../')

const express = require('express')
const app = express()
const http = require('http').Server(app)
const io = require('socket.io')(http)

var clients = {}

app.get('/', (req, res) => res.sendFile(__dirname + '/client.html'))
app.use(express.static('images'))

io.on('connection', (socket) => {
	// Add the client to the active clients object
	clients[socket.id] = socket
	socket.on('disconnect', () => {
		// Remove the client from the active clients object
		delete clients[socket.id]
	})
})

http.listen(3000, () => {
	console.log('Server is running on port 3000.')
	console.log('To see the spectrogram open this in your browser:', 'http://localhost:3000')
})

// Create a new audio engine
var engine = coreAudio.createNewAudioEngine()
engine.setOptions({
	inputChannels: 1,
	outputChannels: 1,
	sampleFormat: 'sampleFormatFloat32',
	interleaved: false,
	fftWindowSize: 1024,
	fftOverlapSize: 0.0,
	fftWindowFunction: 'Blackman'
})

/*function processAudio(samples) {
	return samples;
}
engine.addAudioCallback(processAudio)*/

engine.setFFTCallback((channel, samples) => {
	// Send the fft results to all connected clients.
	for (var id in clients) {
		clients[id].emit('fft', channel, samples)
	}
})