#!/usr/bin/python

"""
For more information on this program, please type
python my_filename.py -h
or
./my_filename.py -h
"""

from subprocess import call,check_output
import re
import inspect
import ast
from tempfile import mkstemp
from shutil import move
from os import remove, close
import json
import vtk
import numpy
import argparse

import pdb

def printStack(startFromTopLevel=0):
    stack = inspect.stack()
    for i in range(len(stack)-1,startFromTopLevel,-1) :
        callerframerecord = stack[i]
        frame = callerframerecord[0]
        info = inspect.getframeinfo(frame)
        indent=" "*(len(stack)-1-i)
        print indent+info.filename+':'+str(info.lineno)+' : '+info.function

def addJsonSubnode(json_str,nodeNameParent,nodeName,nodeValue) :
    endPosNodeNameAbove = 0
    if nodeNameParent != '':
            endPosNodeNameAbove = re.search(nodeNameParent,json_str).end();
    posStartNode = re.search('\{',json_str[endPosNodeNameAbove:]).end() + endPosNodeNameAbove;
    posEndNod = re.search('\}',json_str[posStartNode:]).start() ;

    newJson_str = json_str[:posStartNode]+'\n\"'+nodeName+'\" : \n'+nodeValue
    if posEndNod > 2 :
            newJson_str += ','
    newJson_str += '\n'+ json_str[(posStartNode+1):]
    return newJson_str;

def filesInStringToStringList(stringFiles,inputEnding = '') :
    currentFiles = []
    while len(stringFiles) :
            posEndOfLine = re.search('\n',stringFiles).start()
            currentFile = stringFiles[:posEndOfLine]
            currentFile_re = re.search(inputEnding, currentFile)
            if (len(inputEnding) == 0) or (currentFile_re is not None ):
                    currentFiles.append(currentFile)
            stringFiles = stringFiles[posEndOfLine+1:]
    return currentFiles

def lsToStringList(location = '.',inputEnding = '') :
        filesInFolder = check_output(['ls',location])
        stringList = filesInStringToStringList(filesInFolder,inputEnding)
        if location != '.' :
            for id in range(len(stringList)) :
                stringList[id] = location + stringList[id]
        return stringList

def findToStringList(filename,location = '.') :
        filesInFolder = check_output(['find',location,'-name',filename])
        stringList = filesInStringToStringList(filesInFolder)
        return stringList

def replaceStringInFile(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    for line in old_file:
        new_file.write(line.replace(pattern, subst))
    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

def replaceStringsInFile(file_path,patterns,substitutes) :
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    for line in old_file:
        newLine = line
        for i in range(len(patterns)) :
            pattern = patterns[i]
            subst = substitutes[i]
            newLine = newLine.replace(pattern, subst)
        new_file.write(newLine) 
    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

def getPointsAndConnect_vtkUnstructuredGrid(filename) :
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    data = reader.GetOutput()
    # get the points and the point ids in every cell
    pts_vtk = data.GetPoints()
    pts = numpy.matrix([pts_vtk.GetPoint(i) for i in range(pts_vtk.GetNumberOfPoints())])
    cells = []
    cellTypes = []
    for i in range(data.GetNumberOfCells()):
        currentCell = data.GetCell(i)
        cells.append(numpy.array([currentCell.GetPointId(i) for i in range(currentCell.GetNumberOfPoints())]))
        cellTypes.append(currentCell.GetCellType())
    return pts, cells, numpy.array(cellTypes);

def chopStringAtChar(stringToChop,char,useContentBeforeChar=0) :
    choppedString = stringToChop
    choppedString_re = re.search(char,stringToChop[::-1])
    if choppedString_re is not None :
        posChar = choppedString_re.start()
        choppedString = stringToChop[len(stringToChop)-posChar:]
        if useContentBeforeChar :
            choppedString = stringToChop[:-choppedString_re.end()]
    return choppedString;

def getFilenameWithoutPathAndFilenameEndings(filename) :
    filenameWithoutPath = chopStringAtChar(filename,'/')
    filenameWithoutPathAndFilenameEndings = chopStringAtChar(filenameWithoutPath,'\.',useContentBeforeChar=1)
    return filenameWithoutPathAndFilenameEndings;


def readVariable(filename,variableName,useLineAfterVariable=1,onlyFindOneFit=0) :
	var = []
	f = open(filename,'r')
	for line in f :
		lineContains_variable_re = re.search(variableName,line)
		if  lineContains_variable_re is not None :
			lineAfterVariable = line[lineContains_variable_re.end():]
			curValue_re = re.search('[-+]?[0-9]*\.?[0-9]+',line)
			if useLineAfterVariable :
				curValue_re = re.search('[-+]?[0-9]*\.?[0-9]+',lineAfterVariable)
			if curValue_re is not None and not onlyFindOneFit :
				curValue = ast.literal_eval(curValue_re.group(0))
				var.append(curValue)
			elif curValue_re is not None and onlyFindOneFit :
				var = ast.literal_eval(curValue_re.group(0))
	f.close()
	return var ;

def sortListOfStringsByNumberInString(stringList) :
    numbersInString = []
    for curString in stringList :
        curString_re = re.search('\d{1,}',curString[::-1])	# search from the back, since we want to check the filename not the folders!
        if curString_re is None :
            printStack()
            print "ERROR!!!!!!!!!!!!!!!!!!!!!"
        else :
           numbersInString.append(ast.literal_eval(curString_re.group(0)[::-1]))
    sortedIndexes = [i[0] for i in sorted(enumerate(numbersInString), key=lambda x:x[1])]
    sortedList = [stringList[i] for i in sortedIndexes]
    return sortedList;

def createScene(filenameSceneTemplate,filenameScene,patterns,substitutes) :
    call(['cp',filenameSceneTemplate,filenameScene])
    replaceStringsInFile(filenameScene,patterns,substitutes)
    
def RMSErrorBetweenMeshesWithSameConnect(filename1,filename2) :
    # get the meshes
    pts1, connect1, cellTypes1 = getPointsAndConnect_vtkUnstructuredGrid(filename1)
    pts2, connect2, cellTypes2 = getPointsAndConnect_vtkUnstructuredGrid(filename2)
    # check whether they have the same connect
    if len(connect1) != len(connect2) :
        print "WARNING - the meshes do not have the same connectivity"
        print len(connect1),len(connect2)
        return 0.0
    if len(pts1) != len(pts2) :
        print "WARNING - the meshes do not have the same number of points"
        return 0.0
    # compare the position of the points
    nPoints = len(pts1)
    pts1Vector = numpy.hstack(pts1)
    pts2Vector = numpy.hstack(pts2)
    return numpy.sqrt(numpy.mean(numpy.square(pts1Vector-pts2Vector)))

def parseInput() :
    # for more information on the parser go to
    # https://docs.python.org/2/library/argparse.html#module-argparse
    parser = argparse.ArgumentParser(
        description='Script containing useful functions',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,# ArgumentDefaultsHelpFormatter
        epilog='''Script containing useful functions for the simplified usage of the c++ code and for the evaluation of the convergence. Author: Christoph PAULUS, christophjpaulus@gmail.com''')
    parser.add_argument('-l', dest='listFunctions', action='store_const', default=0, const=1, help='list functions of script')
    args = parser.parse_args()
    return parser,args;

def main() :
    # get the input variables from the parser
    parser,args = parseInput()
    my_filename = parser.prog

    # show the functions if demanded
    if args.listFunctions :
        print "\n\n________________________________________________________________________________"
        print "functions of "+my_filename+"\n\n"
        exec("import " + my_filename[:-3])
        functionList = dir(python_utilities) # HERE WE SHOULD ADAPT IT AS WELL, SUCH THAT THE NAME OF THE FILE CAN BE CHANGED AND IT STILL WORKS
        for curFct in functionList :
            print curFct

if __name__ == "__main__" :
    main()