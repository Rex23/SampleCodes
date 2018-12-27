# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.14-1 replay file
# Internal Version: 2014_06_04-18.11.02 134264
# Run by xiang on Mon Oct 22 16:01:48 2018
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *

def abaqus_gen(seed_size):
    session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=294.501281738281, 
        height=166.658569335938)
    session.viewports['Viewport: 1'].makeCurrent()
    session.viewports['Viewport: 1'].maximize()
    from caeModules import *
    from driverUtils import executeOnCaeStartup
    executeOnCaeStartup()
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=ON)
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(2.0, 2.0))
    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-1']
    p.BaseSolidExtrude(sketch=s, depth=2.0)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.38919, 
        farPlane=9.38376, width=5.47521, height=2.52895, cameraPosition=(5.06331, 
        -0.11133, 7.06797), cameraUpVector=(-0.189834, 0.975435, -0.111759), 
        cameraTarget=(1.04196, 0.916084, 1.04196))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.23348, 
        farPlane=9.60382, width=5.31701, height=2.45588, cameraPosition=(4.37142, 
        -3.15971, 6.13531), cameraUpVector=(-0.165474, 0.918619, 0.358828), 
        cameraTarget=(1.03546, 0.887444, 1.0332))
    p = mdb.models['Model-1'].parts['Part-1']
    f = p.faces
    p.DatumPlaneByOffset(plane=f[4], flip=SIDE2, offset=0.4)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.33512, 
        farPlane=9.62891, width=5.42028, height=2.50358, cameraPosition=(2.72059, 
        -4.98849, -3.14231), cameraUpVector=(0.353513, -0.121716, 0.927477), 
        cameraTarget=(1.01286, 0.862406, 0.906177))
    p = mdb.models['Model-1'].parts['Part-1']
    f1 = p.faces
    p.DatumPlaneByOffset(plane=f1[5], flip=SIDE2, offset=0.4)
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    d1 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=d1[2], cells=pickedCells)
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    pickedCells = c.getSequenceFromMask(mask=('[#1 ]', ), )
    d2 = p.datums
    p.PartitionCellByDatumPlane(datumPlane=d2[3], cells=pickedCells)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.56849, 
        farPlane=9.42069, width=5.65739, height=2.61309, cameraPosition=(-1.02335, 
        -6.17215, 0.202649), cameraUpVector=(-0.0700509, 0.260754, 0.96286), 
        cameraTarget=(0.930327, 0.836313, 0.979914))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.45817, 
        farPlane=9.52608, width=5.54531, height=2.56132, cameraPosition=(-0.31057, 
        -5.95421, 3.46045), cameraUpVector=(0.100623, 0.616468, 0.780924), 
        cameraTarget=(0.94721, 0.841475, 1.05708))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.42645, 
        farPlane=9.56508, width=5.51308, height=2.54643, cameraPosition=(-1.48474, 
        -5.89498, 2.57237), cameraUpVector=(0.159833, 0.49911, 0.851671), 
        cameraTarget=(0.919775, 0.842859, 1.03633))
    session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)
    session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
        referenceRepresentation=OFF)
    p = mdb.models['Model-1'].parts['Part-1']
    p.seedPart(size=seed_size, deviationFactor=seed_size, minSizeFactor=seed_size) #Xiang: seed size
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.12952, 
        farPlane=9.86201, width=9.40242, height=4.35729, viewOffsetX=1.46741, 
        viewOffsetY=-0.680017)
    elemType1 = mesh.ElemType(elemCode=C3D8RH, elemLibrary=STANDARD, 
        kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#2 ]', ), )
    pickedRegions =(cells, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
        elemType3))
    p = mdb.models['Model-1'].parts['Part-1']
    p.generateMesh()
    #: 
    #: Element 7701: Linear hexahedron, type C3D8R
    #:    Nodal connectivity: 2698, 9091, 9110, 2699, 2695, 9034, 9053, 2696
    #: 
    #: Element 1240: Linear hexahedron, type C3D8R
    #:    Nodal connectivity: 3998, 605, 606, 4017, 3941, 602, 603, 3960
    #: 
    #: Element 5321: Linear hexahedron, type C3D8RH
    #:    Nodal connectivity: 1927, 7438, 7457, 1928, 1916, 7229, 7248, 1917
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.25603, 
        farPlane=9.7355, width=7.07066, height=3.27671, viewOffsetX=1.2557, 
        viewOffsetY=-0.166592)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=6.5788, 
        farPlane=10.7986, width=8.85012, height=4.10135, cameraPosition=(3.76526, 
        -7.24544, 1.36227), cameraUpVector=(-0.194094, 0.338803, 0.920619), 
        cameraTarget=(1.046, -0.465354, 0.943132), viewOffsetX=1.57172, 
        viewOffsetY=-0.208517)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=6.34693, 
        farPlane=11.2801, width=8.53819, height=3.95679, cameraPosition=(4.65722, 
        -6.71285, 3.24768), cameraUpVector=(-0.195882, 0.521114, 0.830705), 
        cameraTarget=(1.26464, -0.477612, 1.47236), viewOffsetX=1.51632, 
        viewOffsetY=-0.201168)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=6.34825, 
        farPlane=11.186, width=8.53996, height=3.95761, cameraPosition=(4.39278, 
        -6.82883, 3.07406), cameraUpVector=(-0.213771, 0.501278, 0.838464), 
        cameraTarget=(1.19309, -0.473529, 1.3677), viewOffsetX=1.51663, 
        viewOffsetY=-0.20121)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=6.37265, 
        farPlane=11.4327, width=8.57279, height=3.97282, cameraPosition=(5.18011, 
        -6.20571, 4.18012), cameraUpVector=(-0.309084, 0.570225, 0.761124), 
        cameraTarget=(1.39853, -0.506918, 1.57954), viewOffsetX=1.52246, 
        viewOffsetY=-0.201983)
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    p = mdb.models['Model-1'].parts['Part-1']
    a.Instance(name='Part-1-1', part=p, dependent=ON)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        adaptiveMeshConstraints=ON)
    mdb.models['Model-1'].StaticStep(name='Step-1', previous='Initial')
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
    mdb.models['Model-1'].steps['Step-1'].setValues(timePeriod=1.0, maxNumInc=1000, minInc=1.0)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.61264, 
        farPlane=10.1796, width=5.18423, height=2.65817, cameraPosition=(5.42492, 
        -5.09859, 3.36137), cameraUpVector=(-0.157365, 0.616011, 0.771859))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.66988, 
        farPlane=10.1259, width=5.2371, height=2.68528, cameraPosition=(4.23869, 
        -5.67106, 3.71751), cameraUpVector=(-0.160216, 0.627072, 0.762307), 
        cameraTarget=(1.03404, 0.908058, 1.03953))
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
        predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
    a = mdb.models['Model-1'].rootAssembly
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#2000 ]', ), )
    region = a.Set(vertices=verts1, name='Set-1')
    mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Step-1', 
        region=region, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-1'].rootAssembly
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#4000 ]', ), )
    region = a.Set(vertices=verts1, name='Set-2')
    mdb.models['Model-1'].DisplacementBC(name='BC-2', createStepName='Step-1', 
        region=region, u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-1'].rootAssembly
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#200 ]', ), )
    region = a.Set(vertices=verts1, name='Set-3')
    mdb.models['Model-1'].DisplacementBC(name='BC-3', createStepName='Step-1', 
        region=region, u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    a = mdb.models['Model-1'].rootAssembly
    v1 = a.instances['Part-1-1'].vertices
    verts1 = v1.getSequenceFromMask(mask=('[#10 ]', ), )
    region = a.Set(vertices=verts1, name='Set-4')
    mdb.models['Model-1'].DisplacementBC(name='BC-4', createStepName='Step-1', 
        region=region, u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.6084, 
        farPlane=10.1817, width=5.18032, height=2.65616, cameraPosition=(3.97928, 
        -5.04982, -3.10554), cameraUpVector=(0.0731364, -0.193325, 0.978405), 
        cameraTarget=(1.03058, 0.916333, 0.948642))
    mdb.models['Model-1'].TabularAmplitude(name='Amp-1', timeSpan=STEP, 
        smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1e-20, 1.0), (1.0, 1.0)))
    a = mdb.models['Model-1'].rootAssembly
    s1 = a.instances['Part-1-1'].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#c000 ]', ), )
    region = a.Surface(side1Faces=side1Faces1, name='Surf-1')
    mdb.models['Model-1'].Pressure(name='Load-1', createStepName='Step-1', 
        region=region, distributionType=UNIFORM, field='', magnitude=-1000000.0, 
        amplitude='Amp-1')
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.60794, 
        farPlane=10.1846, width=5.17989, height=2.65594, cameraPosition=(4.61703, 
        -5.19809, 4.2942), cameraUpVector=(-0.501417, 0.50495, 0.702572), 
        cameraTarget=(1.03885, 0.914411, 1.04457))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.83078, 
        farPlane=9.95801, width=5.38572, height=2.76148, cameraPosition=(5.82833, 
        -5.24412, 1.14066), cameraUpVector=(-0.182822, 0.300993, 0.935938), 
        cameraTarget=(1.05474, 0.913807, 1.00321))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.79164, 
        farPlane=9.99864, width=5.34957, height=2.74294, cameraPosition=(5.23577, 
        -5.64333, 0.492364), cameraUpVector=(-0.145533, 0.22961, 0.962341), 
        cameraTarget=(1.04711, 0.908664, 0.994858))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.75282, 
        farPlane=10.0287, width=5.31372, height=2.72456, cameraPosition=(7.91752, 
        -2.59847, 2.20949), cameraUpVector=(-0.470511, 0.119482, 0.874267), 
        cameraTarget=(1.08191, 0.948175, 1.01714))
    session.viewports['Viewport: 1'].view.setValues(nearPlane=5.83926, 
        farPlane=9.94034, width=5.39356, height=2.7655, cameraPosition=(8.29367, 
        -1.85364, 1.95298), cameraUpVector=(-0.422925, 0.138868, 0.895461), 
        cameraTarget=(1.08659, 0.957433, 1.01395))
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=OFF, bcs=OFF, 
        predefinedFields=OFF, connectors=OFF)
    mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
        numGPUs=0)
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
        engineeringFeatures=ON, mesh=OFF)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=OFF)
    p1 = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#5 ]', ), )
    p.Set(cells=cells, name='Abaqus')
    #: The set 'Abaqus' has been created (2 cells).
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#2 ]', ), )
    p.Set(cells=cells, name='XFEM')
    #: The set 'XFEM' has been created (1 cell).
    mdb.models['Model-1'].Material(name='Material-XFEM')
    mdb.models['Model-1'].materials['Material-XFEM'].Elastic(table=((1000000000.0, 
        0.3), ))
    mdb.models['Model-1'].HomogeneousSolidSection(name='Section_XFEM', 
        material='Material-XFEM', thickness=None)
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    cells = c.getSequenceFromMask(mask=('[#2 ]', ), )
    region = p.Set(cells=cells, name='Set-4')
    p = mdb.models['Model-1'].parts['Part-1']
    p.SectionAssignment(region=region, sectionName='Section_XFEM', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    mdb.models['Model-1'].Material(name='Material-Abaqus')
    mdb.models['Model-1'].materials['Material-Abaqus'].Elastic(table=((
        1000000000.0, 0.3), ))
    mdb.models['Model-1'].HomogeneousSolidSection(name='Section_Abaqus', 
        material='Material-Abaqus', thickness=None)
    p = mdb.models['Model-1'].parts['Part-1']
    region = p.sets['Abaqus']
    p = mdb.models['Model-1'].parts['Part-1']
    p.SectionAssignment(region=region, sectionName='Section_Abaqus', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    p1 = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p1 = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p1 = mdb.models['Model-1'].parts['Part-1']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    mdb.jobs['Job-1'].writeInput(consistencyChecking=OFF)
    #: The job input file has been written to "Job-1.inp".
    mdb.saveAs(pathName='./Test')
    #: The model database has been saved to "C:\Users\xiang\Documents\Xiang\Papers\Study_SIFs\Test.cae".
    mdb.save()
    #: The model database has been saved to "C:\Users\xiang\Documents\Xiang\Papers\Study_SIFs\Test.cae".
