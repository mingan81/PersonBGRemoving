#!/usr/bin/python -u

import videoProcessing 
import sys

pathName = sys.argv[1]
pathCommand = sys.argv[2]
pathModel = sys.argv[3]

vd = videoProcessing.videoDecode(pathName,pathCommand,pathModel)
vd.audioProcessing()
vd.videoProcessig()

ve = videoProcessing.videoEncode(vd.videoName)
ve.videoEncoding(17,"webm")
ve.audioEncoding()

vd.deleteFile()
