from __future__ import division
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import wave, struct
import math



def init(fileName):
    # Define global variables
    global frame
    global signal
    global length
    
    # Extract audio from video file
    command = "ffmpeg -y -i " + fileName + " -ab 160k -ac 2 -ar 44100 -vn audio.wav"

    subprocess.call(command, shell=True)

    # Open audio file
    #spf = wave.open('/Users/kevinwang264/Desktop/test.wav','r')
    spf = wave.open('audio.wav','r')
    frame = spf.getframerate()*2

    # Extract Raw Audio data from Wav File
    sig = spf.readframes(-1)
    signal = np.fromstring(sig, 'Int16')

    # length of audio in frames
    length = len(signal) 


def getFirstMaxima():
    # Find first maxima corresponding to start
    max = 0
    for i in range(0,int(length/2)):
        if signal[i] > max:
            max = signal[i]

    t = 0
    for i in range(0,int(length/2)):
        if signal[i] == max:
            t = i

    maxima = t/frame
    return maxima


def getSecondMaxima():
    # Find second maxima corresponding to stop/end
    max = 0
    for i in range(int(length/2), length):
        if signal[i] > max:
            max = signal[i]

    t = 0
    for i in range(int(length/2), length):
        if signal[i] == max:
            t = i

    maxima = t/frame
    return maxima


def getVideoLength(fileName):

    init(fileName)

    firstMaxima = getFirstMaxima()
    print firstMaxima
    secondMaxima = getSecondMaxima()
    print secondMaxima
    videoLength = secondMaxima - firstMaxima

    return videoLength;





