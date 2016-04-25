import subprocess

subprocess.call("""
	sddsanaylzebeam input.bun input.anb""",shell=True)

subprocess.call("""
	sddsanalyzebeam output.bun output.anb""",shell=True)
emitsIn = subprocess.Popen("""
sdds2stream input.anb -col=enx,eny -delimiter=" "
""", shell = True, stdout = subprocess.PIPE)

emitsOut = subprocess.Popen("""
sdds2stream output.anb -col=enx,eny -delimiter=" "
""", shell = True, stdout = subprocess.PIPE)

inVals = emitsIn.stdout.read()
enx0 = float(inVals.split()[0])
eny0 = float(inVals.split()[1])

outVals = emitsOut.stdout.read()
enxf = float(outVals.split()[0])
enyf = float(outVals.split()[1])


print "enx0 = %s enxf = %s f/0 = %s\neny0 = %s enyf = %s f/0 = %s\n" % (enx0,enxf,enxf/enx0,eny0,enyf,enyf/eny0)