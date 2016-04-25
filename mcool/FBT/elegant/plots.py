import subprocess
#-scale=-0.015,0.015,-0.015,0.015
subprocess.call("""
	sddsplot input.bun -title='Input X-Y (Flat Beam)' \
	-col=x,y,xp,yp -scale=-0.022,0.022,-0.022,0.022 \
	-arrowSettings=autoscale,cartesianData""",shell=True)

subprocess.call("""
	sddsplot output.bun -title='Output X-Y (FBT-Solenoid-IFBT)' \
	-col=x,y,xp,yp -scale=-0.022,0.022,-0.022,0.022 \
	-arrowSettings=autoscale,cartesianData""",shell=True)


