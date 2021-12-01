#!/usr/bin/python

"""
For more information on this program, please type
python my_filename.py -h
or
./my_filename.py -h
"""
from os import chdir
from subprocess import call,check_output
import argparse
from python_utilities import findToStringList, replaceStringInFile, getFilenameWithoutPathAndFilenameEndings
import pdb

def parseInput() :
    # for more information on the parser go to
    # https://docs.python.org/2/library/argparse.html#module-argparse
    parser = argparse.ArgumentParser(
        description='Script to simplify the usage of the xfem implementation in c++',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,# ArgumentDefaultsHelpFormatter
        epilog='''This script has been tested on Ubuntu 14.04. Author: Christoph PAULUS, christophjpaulus@gmail.com''')
    parser.add_argument('-l', dest='installLibraries', action='store_const', default=0, const=1, help='install necessary libraries')
    parser.add_argument('--cTINL', metavar='C', nargs='*', type=str, default=['sudo','apt-get','install'], help='command to install new libraries') 
    parser.add_argument('-c', dest='compile', action='store_const', default=0, const=1, help='compile code')
    parser.add_argument('--cf', metavar='CF', nargs='?', default='build', help='folder for the compilation')
    parser.add_argument('-p', dest='setPath', action='store_const', default=0, const=1, help='set correct paths in the xml files')
    parser.add_argument('-o', dest='createOutputFolders', action='store_const', default=0, const=1, help='create output folders for each xml scene')
    parser.add_argument('-e', dest='lanceExample', action='store_const', default=0, const=1, help='lance example')
    parser.add_argument('--example', metavar='E', nargs='+', default='scenes/Liver.xml', help='example(s) to lance')
    parser.add_argument('-a', dest='lanceEvaluation', action='store_const', default=0, const=1, help='lance evaluation/analysis of the implementation')
    args = parser.parse_args()
    return parser,args;

def main() :
    # get the input variables from the parser
    parser,args = parseInput()
    compilationFolder = args.cf
    example = args.example

    # get all the xml files in the current repository
    xmlFilenames = findToStringList('*.xml')

    # installation of the libraries
    if args.installLibraries :
        print "\n\n________________________________________________________________________________"
        print "install necessary libraries\n\n"
        libraryList = ['libeigen3-dev','libvtk5-dev','libcgal-dev','python-matplotlib','libtinyxml-dev','paraview','cmake','make']
        for curLib in libraryList :
            args.cTINL.append(curLib)
            call(args.cTINL)
            args.cTINL.remove(curLib)

    # compilation of the xfem code
    if args.compile :
        print "\n\n________________________________________________________________________________"
        print "compilation of the xfem code in the folder "+compilationFolder+"\n\n"
        call(['mkdir',compilationFolder])
        currentPath = check_output(['pwd'])[:-1]
        chdir(currentPath+'/'+ compilationFolder)
        call(['cmake','..'])
        call(['make'])
        chdir(currentPath)

    # set the correct path in the existing xml files
    if args.setPath :
        print "\n\n________________________________________________________________________________"
        print "set correct paths in the xml files\n\n"
        currentPath = check_output(['pwd'])[:-1]
        for curXml in xmlFilenames :
            print curXml
            replaceStringInFile(curXml,'PATH',currentPath)

    # introduction of the necessary output folders
    if args.createOutputFolders :
        print "\n\n________________________________________________________________________________"
        print "create output folders for each xml scene\n\n"
        for curXml in xmlFilenames :
            print curXml
            sceneName = getFilenameWithoutPathAndFilenameEndings(curXml)
            call(['mkdir','-p','output/'+sceneName+'/displacedObject'])
            call(['mkdir','-p','output/'+sceneName+'/stiffness'])

    # lance the example Liver.xml
    if args.lanceExample :
        print "\n\n________________________________________________________________________________"
        print "lance the example "+example+"\n\n"
        call([compilationFolder+'/XFEMExec',example])
        call(['paraview','--data=output/Liver/displacedObject/XFEM_Visualization.vtk'])

    # lance evaluation/analysis of the implementation
    if args.lanceEvaluation :
        print "\n\n________________________________________________________________________________"
        print "lance evaluation/analysis of the implementation\n\n"
        call(['python','evaluateConvergence.py','FEM.json','XFEM.json','-c',compilationFolder])

if __name__ == "__main__" :
    main()