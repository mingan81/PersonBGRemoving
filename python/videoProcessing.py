import sys
import subprocess


class videoDecode:
	def __init__(self,pathName,pathCommand,pahtModel):
		pathSplit = pathName.split("/")
		self.pathName = pathName
		self.pathCommand = pathCommand
		self.pahtModel = pahtModel
		self.videoName = pathSplit[len(pathSplit) - 1].split(".")[0]
		print self.videoName
		#command = "ffmpeg -i " + self.pathName + " 2>&1|sed -n 's/.*, \(.*\) fp.*/\1/p'"
		#self.fps = subprocess.check_output(command, shell=True)
		#print type(self.fps)
		#print float(self.fps)
		#command = "ffprobe -i " + self.pathName + " 2>&1 | grep Stream | grep -oP ', \K[0-9]+x[0-9]+'"
		#self.resolution = subprocess.check_output(command, shell=True)
		#self.w = self.resolution.split("x")[0]
		#self.h = self.resolution.split("x")[1]
		
	def audioProcessing(self):
		command = "ffmpeg -y -i " + self.pathName + " -ab 160k -ac 2 -ar 44100 -vn " + self.videoName + ".wav"
		subprocess.call(command, shell=True)
	
	def videoProcessig(self):
		command = "mkdir " + self.videoName
		subprocess.call(command, shell=True)
		
		print self.pathName
		print self.videoName
		command = "ffmpeg -i " + self.pathName +" " + self.videoName + "/" + self.videoName  + "%04d.png"
		subprocess.call(command, shell=True)
		
		command = self.pathCommand + " " +  self.videoName + "/ " +  self.videoName + " " + self.pahtModel
		subprocess.call(command, shell=True)

	def deleteFile(self):
		command = "rm " + self.videoName + ".wav"
		subprocess.call(command, shell=True)
		
		command = "rm -r " + self.videoName
		subprocess.call(command, shell=True)
		
		

class videoEncode:
	def __init__(self,videoName):
		#~ pathSplit = pathName.split("/")
		#~ self.videoName = pathSplit[len(pathSplit) - 2]
		self.videoName = videoName
	
	def videoEncoding(self,fps,videoFormat):
		self.fps = fps
		self.videoFormat = videoFormat
		command = "ffmpeg -r " + str(self.fps) + " -i " + self.videoName + "/" + self.videoName +"_out%04d.png " + self.videoName +"_outWoAudio." + self.videoFormat
		subprocess.call(command, shell=True)
	
	def audioEncoding(self):
		command = "ffmpeg " + "-i " + self.videoName + "_outWoAudio." + self.videoFormat + " -i " + self.videoName + ".wav " + self.videoName +"_out." + self.videoFormat
		subprocess.call(command, shell=True)
	
		
