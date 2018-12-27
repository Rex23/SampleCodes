# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.14-1 replay file
# Internal Version: 2014_06_04-18.11.02 134264
# Run by xiang on Tue Oct 23 22:00:25 2018
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=251.29817199707, 
    height=166.658569335938)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup

def generate_crack(theta):
	executeOnCaeStartup()
	session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
		referenceRepresentation=ON)
	s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
		sheetSize=200.0)
	g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
	s.setPrimaryObject(option=STANDALONE)
	s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.15, 0.0))
	p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
		type=DEFORMABLE_BODY)
	p = mdb.models['Model-1'].parts['Part-1']
	p.BaseShell(sketch=s)
	s.unsetPrimaryObject()
	p = mdb.models['Model-1'].parts['Part-1']
	session.viewports['Viewport: 1'].setValues(displayedObject=p)
	del mdb.models['Model-1'].sketches['__profile__']
	session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
	session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
		meshTechnique=ON)
	session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
		referenceRepresentation=OFF)
	p = mdb.models['Model-1'].parts['Part-1']
	p.seedPart(size=0.02, deviationFactor=0.1, minSizeFactor=0.1)
	p = mdb.models['Model-1'].parts['Part-1']
	p.generateMesh()
	session.viewports['Viewport: 1'].view.setValues(nearPlane=0.778122, 
		farPlane=0.918934, width=0.764818, height=0.393598, viewOffsetX=0.0594054, 
		viewOffsetY=-0.0115227)
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
	p.deleteMesh(regions=pickedRegions)
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
	p.setMeshControls(regions=pickedRegions, elemShape=TRI)
	p = mdb.models['Model-1'].parts['Part-1']
	p.generateMesh()
	session.viewports['Viewport: 1'].view.setValues(nearPlane=0.774, 
		farPlane=0.923056, width=0.809327, height=0.416503, viewOffsetX=0.0729637, 
		viewOffsetY=-0.0143696)
	a = mdb.models['Model-1'].rootAssembly
	session.viewports['Viewport: 1'].setValues(displayedObject=a)
	session.viewports['Viewport: 1'].assemblyDisplay.setValues(
		optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
	a = mdb.models['Model-1'].rootAssembly
	a.DatumCsysByDefault(CARTESIAN)
	p = mdb.models['Model-1'].parts['Part-1']
	a.Instance(name='Part-1-1', part=p, dependent=ON)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=0.7603, 
		farPlane=0.936756, width=0.960684, height=0.492582, viewOffsetX=0.0476149, 
		viewOffsetY=0.0566029)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=0.54053, 
		farPlane=1.03279, width=0.682992, height=0.350198, cameraPosition=(
		0.0055824, -0.787669, -0.0249881), cameraUpVector=(0.260959, 0.0313231, 
		0.964841), cameraTarget=(-0.0146282, 0.0603319, -0.0470517), 
		viewOffsetX=0.0338515, viewOffsetY=0.0402414)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=0.598516, 
		farPlane=0.979142, width=0.75626, height=0.387766, cameraPosition=(
		0.0842667, -0.630675, 0.471793), cameraUpVector=(-0.258472, 0.62848, 
		0.733625), cameraTarget=(0.0121446, -0.00132998, -0.0927631), 
		viewOffsetX=0.0374829, viewOffsetY=0.0445583)
	a = mdb.models['Model-1'].rootAssembly
	a.DatumPointByCoordinate(coords=(1.0, 1.0, 1.0))
	session.viewports['Viewport: 1'].view.setValues(nearPlane=1.61175, 
		farPlane=3.71504, width=2.64841, height=1.35795, viewOffsetX=0.474637, 
		viewOffsetY=0.35374)
	a1 = mdb.models['Model-1'].rootAssembly
	a1.translate(instanceList=('Part-1-1', ), vector=(1.0, 1.0, 1.0))
	#: The instance Part-1-1 was translated by 1., 1., 1. with respect to the assembly coordinate system
	session.viewports['Viewport: 1'].view.setValues(nearPlane=1.52388, 
		farPlane=3.99993, width=4.64898, height=2.38372, viewOffsetX=1.26953, 
		viewOffsetY=0.349317)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=1.90233, 
		farPlane=4.66285, width=5.80354, height=2.97571, cameraPosition=(1.61314, 
		-2.56178, 1.68152), cameraUpVector=(-0.382872, 0.457139, 0.802766), 
		cameraTarget=(0.418795, -0.671365, 0.0353838), viewOffsetX=1.58481, 
		viewOffsetY=0.436069)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=2.2275, 
		farPlane=5.20247, width=6.79554, height=3.48435, cameraPosition=(2.63159, 
		-2.70834, 0.93819), cameraUpVector=(-0.166742, 0.250763, 0.95358), 
		cameraTarget=(0.707107, -0.881118, 0.121173), viewOffsetX=1.8557, 
		viewOffsetY=0.510606)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=1.75191, 
		farPlane=4.92159, width=5.34463, height=2.74041, cameraPosition=(1.81713, 
		-2.73711, 0.878547), cameraUpVector=(-0.210758, 0.266733, 0.940444), 
		cameraTarget=(0.285064, -0.622155, -0.0646495), viewOffsetX=1.45949, 
		viewOffsetY=0.401587)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=1.93412, 
		farPlane=4.84631, width=5.9005, height=3.02543, cameraPosition=(1.88819, 
		-2.70719, 1.15949), cameraUpVector=(-0.213735, 0.323084, 0.921919), 
		cameraTarget=(0.321044, -0.681776, 0.0863623), viewOffsetX=1.61129, 
		viewOffsetY=0.443355)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=2.27049, 
		farPlane=4.50994, width=1.77557, height=0.910407, viewOffsetX=1.36701, 
		viewOffsetY=0.8291)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=2.30603, 
		farPlane=4.52222, width=1.80337, height=0.924661, cameraPosition=(1.25783, 
		-2.56805, 2.14858), cameraUpVector=(-0.200483, 0.561696, 0.802686), 
		cameraTarget=(0.0569273, -0.667665, 0.518794), viewOffsetX=1.38841, 
		viewOffsetY=0.842081)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=2.3337, 
		farPlane=4.57829, width=1.82501, height=0.935754, cameraPosition=(1.57662, 
		-2.60662, 1.991), cameraUpVector=(-0.201119, 0.516815, 0.832137), 
		cameraTarget=(0.181063, -0.733188, 0.490161), viewOffsetX=1.40507, 
		viewOffsetY=0.852183)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=2.36048, 
		farPlane=4.5515, width=1.44123, height=0.738976, viewOffsetX=1.34173, 
		viewOffsetY=0.855095)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=2.23597, 
		farPlane=4.34873, width=1.3652, height=0.699996, cameraPosition=(0.989445, 
		-2.66775, 1.75854), cameraUpVector=(-0.20532, 0.493506, 0.84516), 
		cameraTarget=(-0.0340327, -0.554781, 0.276074), viewOffsetX=1.27096, 
		viewOffsetY=0.80999)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=2.21219, 
		farPlane=4.31718, width=1.35069, height=0.692552, viewOffsetX=1.32665, 
		viewOffsetY=0.859948)
	a1 = mdb.models['Model-1'].rootAssembly
	a1.rotate(instanceList=('Part-1-1', ), axisPoint=(1.0, 1.0, 1.0), 
		axisDirection=(0.15, 0.0, 0.0), angle=theta)
	#: The instance Part-1-1 was rotated by theta degrees about the axis defined by the point 1., 1., 1. and the vector 150.E-03, 0., 0.
	session.viewports['Viewport: 1'].view.setValues(nearPlane=2.17036, 
		farPlane=4.35901, width=1.80561, height=0.925811, viewOffsetX=1.43557, 
		viewOffsetY=0.871893)
	session.viewports['Viewport: 1'].view.setValues(nearPlane=3.38764, 
		farPlane=5.25722, width=2.81831, height=1.44506, cameraPosition=(5.23305, 
		-0.314876, 0.812927), cameraUpVector=(-0.40101, -0.249389, 0.881474), 
		cameraTarget=(2.71319, -0.378972, -0.351586), viewOffsetX=2.24073, 
		viewOffsetY=1.3609)
	mdb.Job(name='Cracks_Ori', model='Model-1', description='', type=ANALYSIS, 
		atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
		memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
		explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
		modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
		scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, 
		numDomains=1, activateLoadBalancing=False, multiprocessingMode=DEFAULT, 
		numCpus=1)
	mdb.jobs['Cracks_Ori'].writeInput(consistencyChecking=OFF)
	mdb.saveAs(
		pathName='./Cracks.cae')
	#: The model database has been saved to "C:\Users\xiang\Documents\Xiang\Papers\Study_SIFs\Crack\Crack.cae".
	mdb.save()
	#: The model database has been saved to "C:\Users\xiang\Documents\Xiang\Papers\Study_SIFs\Crack\Crack.cae".
