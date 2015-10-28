from __future__ import division
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import wave, struct
import math

#Extract audio from video file
command = "ffmpeg -i /home/rajkumar/rpv.mkv -ab 160k -ac 2 -ar 44100 -vn audio.wav"

subprocess.call(command, shell=True)

#open text file to write time values
f = open('time.txt','wb')

#open audio file
spf = wave.open('audio.wav','r')
frame = spf.getframerate()*2

#Extract Raw Audio data from Wav File
signal = spf.readframes(-1)
signal = np.fromstring(signal, 'Int16')

#length of audio in frames
length = len(signal) 

#find first maxima corresponding to start
max = 0
for i in range(0,int(length/2)):
	if signal[i] > max:
		max = signal[i]

t = 0
for i in range(0,int(length/2)):
	if signal[i] == max:
		t = i

x = t/frame

s = str(x)

f.write(s + '\n')

#find second maxima corresponding to stop/end
max = 0
for i in range(int(length/2), length):
	if signal[i] > max:
		max = signal[i]

t = 0
for i in range(int(length/2), length):
	if signal[i] == max:
		t = i

x = 0.0
x = t/frame

s = str(x)
f.write(s)

#Optional if you want to plot the signal 
#Time=np.linspace(0,len(signal)/frame, num=len(signal))

#plt.figure(1)
#plt.title('Signal Wave...')
#plt.plot(Time,signal)
#plt.show()
