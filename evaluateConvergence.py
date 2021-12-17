#!/usr/bin/python

"""
For more information on this program, please type
python my_filename.py -h
or
./my_filename.py -h
"""
import argparse
from subprocess import call, check_output
import json
from python_utilities import createScene,lsToStringList,chopStringAtChar,readVariable,findToStringList,sortListOfStringsByNumberInString,addJsonSubnode,getPointsAndConnect_vtkUnstructuredGrid,RMSErrorBetweenMeshesWithSameConnect
import re
import numpy
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from datetime import datetime

import pdb

def parseInput() :
    # for more information on the parser go to
    # https://docs.python.org/2/library/argparse.html#module-argparse
    parser = argparse.ArgumentParser(
        description='Script for the evaluation of the xfem implementation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,# ArgumentDefaultsHelpFormatter
        epilog='''This script has been tested on Ubuntu 14.04. Author: Christoph PAULUS, christophjpaulus@gmail.com''')
    parser.add_argument('jsonFiles',metavar='J',type=str,nargs='+',help='json files containing the information for each line')
    parser.add_argument('-l', metavar='L', nargs='*', type=int, default=[10,5000], help='limits for the number of degrees of freedom')
    parser.add_argument('-o', nargs='?', default='evaluation', help='filename for output figure')
    parser.add_argument('-c', nargs='?', default='build', help='folder of compiled code')
    parser.add_argument('-p', dest='plotOnly', action='store_const', default=0, const=1, help='plot only, no calculations')
    args = parser.parse_args()
    return parser,args;

def main() :
    # get the input variables from the parser
    parser,args = parseInput()
    minimalNDOF = args.l[0]
    maximalNDOF = args.l[1]
    filenamesOptions = args.jsonFiles
    # preparation of the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel("NDoFs")
    ax.set_ylabel("$\epsilon_{RMS}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    # run over json files containing the information for each line
    for filenameOptions in filenamesOptions :
        # get all information stored in json
        f = open(filenameOptions,'r')
        options_json = json.load(f)
        f.close()
        name = options_json['name']
        calculateDeformations = options_json['calculateDeformations']
        folderObjects = options_json['folderObjects']
        numberOfIntegrationPoints = options_json['numberOfIntegrationPoints']
        filenameCut = options_json['filenameCut']
        sceneTemplate = options_json['sceneTemplate']['filename']
        patterns = options_json['sceneTemplate']['patterns']
        variablesToRead = options_json['variablesToRead']
        filenameReferenceObjectInitialPos = options_json['filenameReferenceObjectInitialPos']
        deformedReferenceMeshName = options_json['deformedReferenceMeshName']
        referenceSolution = options_json['referenceSolution']
        plotErrorLine = options_json['plotErrorLine']
        plotLineColor = options_json['plotLineColor']
        # get other variables
        outputFolder = 'output_'+name+'/'
        currentPath = check_output(['pwd'])[:-1]+'/'
        scenes = 'scenes/'
        pathScenes  = currentPath+scenes
        filenameReferenceObjectInitialPos = currentPath+filenameReferenceObjectInitialPos   
        # output
        print("\n\n"+filenameOptions)
        print("_______________________________"+name+"___________________________________________")

        # calculation of the deformation
        if calculateDeformations and not args.plotOnly :
            print("_______________________________"+name+": Calculation start"+"________________________")
            print(datetime.now())
            # get the executable
            executable = check_output(['find', args.c, '-name', 'XFEMExec'])[:-1]
            # get existent meshes
            pathMeshObject = currentPath+scenes+folderObjects
            objectVTKs = lsToStringList(location = pathMeshObject,inputEnding = '.vtk')
            sortedObjectVTKs = sortListOfStringsByNumberInString(objectVTKs)
            filenameCut = currentPath+filenameCut
            # for each vtk in folderObjects, we create a scene, lance it and analyse its output
            for curObjectVTK in sortedObjectVTKs :
        #        curObjectVTK_re = re.search('Beam_', curObjectVTK)
        #        if curObjectVTK_re is not None :
                    # write new scene file
                    currentFile = chopStringAtChar(curObjectVTK,'/')
                    curNDOF = 3*readVariable(pathMeshObject+currentFile,'POINTS',onlyFindOneFit=1)
                    if minimalNDOF < curNDOF and curNDOF < maximalNDOF :
                        # creation of the scene
                        filenameOutput = currentPath+outputFolder+currentFile[:-4]+'/'
                        substitutes = [pathMeshObject+currentFile,str(numberOfIntegrationPoints),filenameOutput,filenameCut,filenameReferenceObjectInitialPos]
                        currentSceneFilename = pathScenes+currentFile[:-4]+'.xml'
                        createScene(sceneTemplate,currentSceneFilename,patterns,substitutes)
                        # execution of the scene
                        call(['mkdir', '-p', filenameOutput+'displacedObject/'])
                        print(executable+' '+currentSceneFilename)
                        print("nDOF "+str(curNDOF))
                        print(datetime.now())
                        executableOutputFilename = filenameOutput+'executableOutput'
                        f_executableOutput = open(executableOutputFilename,'w') 
                        call([executable,currentSceneFilename],stdout=f_executableOutput)
                        f_executableOutput.close()
                        json_str = "{\n}"
                        # save specific information of execution output in Info.json
                        for curVariableName in variablesToRead :
                            curVariableValue = str(readVariable(executableOutputFilename,curVariableName,onlyFindOneFit=1))
                            json_str = addJsonSubnode(json_str,'',curVariableName,curVariableValue)
                        nodeName = 'deformedHighResolutionMesh'
                        nodeValue = '\"'+filenameOutput+'displacedObject/HighRes.vtk\"'
                        json_str = addJsonSubnode(json_str,'',nodeName,nodeValue)
                        filenameOutput = filenameOutput+'Info.json'
                        f_output = open(filenameOutput,'w')
                        f_output.write(json_str)
                        f_output.close()
            print(datetime.now())
            print("_______________________________"+name+": Calculation end"+"______________________________")

        # comparison of the reference solution with our solution
        if plotErrorLine :
            print("_______________________________"+name+": Plot start"+"_______________________________")
            print(datetime.now())
            # get all the Info.json files that have been written in the outputfolder
            infoFiles = findToStringList('Info.json',outputFolder)
            sortedInfoFiles = sortListOfStringsByNumberInString(infoFiles)
            x = []
            y = []
            for currentFile in sortedInfoFiles :
                # load the content of the Info.json file
                f_InfoJson = open(currentFile,'r')
                InfoJson = json.load(f_InfoJson)
                f_InfoJson.close()
                nDOFs = InfoJson["numberOfDOFs"]
                deformedHighResolutionMesh = InfoJson["deformedHighResolutionMesh"]
                # compare the reference mesh with the solution calculated by our implementation
                print(currentFile)
                print(nDOFs)
                RMS = RMSErrorBetweenMeshesWithSameConnect(currentPath+referenceSolution,deformedHighResolutionMesh)
                print(RMS)
                x.append(nDOFs)
                y.append(RMS)
            # add a new line on the plot
            ax.plot(x,y, '3-'+plotLineColor, label=name)
            plt.legend( loc='upper right')
            print(datetime.now())
            print("_______________________________"+name+": Plot end"+"_________________________________")

    # put the plot into a file
    # plt.show()
    matplotlib.pyplot.savefig(args.o)

if __name__ == "__main__" :
    main()
