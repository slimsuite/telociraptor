#!/usr/bin/python

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
#  
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not, write to 
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       Telociraptor
Description:  Telomere Prediction and Genome Assembly Editing Tool
Version:      0.5.0
Last Edit:    30/04/23
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    Telociraptor is a renamed version of GenomeTweak that also has the capacity for standalone telomere prediction (as
    executed by Diploidocus) along with a combined `chromsyn=T` mode for generating the "telonull" telomere output and
    assembly gaps table used as input by [ChromSyn](https://github.com/slimsuite/chromsyn). The original GenomeTweak
    function is now invoked with `tweak=T`. If running in `tweak` mode, `telomeres` and/or `chromsyn` output will be
    executed on the `seqout=FILE` if generated, which will fixed version of the assembly if `autofix=T`.

    ### ~ Telomere finding [telomeres=T] ~ ###

    Telociraptor performs a regex-based search for Telomeres, based on [FindTelomeres](https://github.com/JanaSperschneider/FindTelomeres).
    By default, this looks for a canonical telomere motif of TTAGGG/CCCTAA, allowing for some variation. Telociraptor
    searches for a forward telomere regex sequence of `C{2,4}T{1,2}A{1,3}` at the 5' end, and a reverse sequence at the
    3' end of `T{1,3}A{1,2}G{2,4}`. These can be set with `telofwd=X` and `telorev=X`. For each sequence, Telociraptor
    trims off any trailing Ns and then searches for telomere-like sequences at sequence ends. For each sequence, the
    presence/absence and length of trimming are reported for the 5' end (tel5 and trim5) and 3' end (tel3 and trim3),
    along with the total percentage telomeric sequence (TelPerc).

    Telomeres are marked if at least 50% (`teloperc=PERC`) of the terminal 50 bp (`telosize=INT`) matches the appropriate
    regex. If either end contains a telomere, the total percentage of the sequence matching either regex is calculated as
    TelPerc. Note that this number neither restricts matches to the termini, not includes sequences within predicted
    telomeres that do not match the regex. By default, only sequences with telomeres are output to the `*.telomeres.tdt`
    output, but switching `telonull=T` will output all sequences. This can be useful if you also need a table of sequence
    lengths, and is the recommended input for [ChromSyn](https://github.com/slimsuite/chromsyn).

    ### ~ GenomeTweak mode [tweak=T] ~ ###

    GenomeTweak is a simple tool designed to help manual curation of genome assembly scaffolding. A genome assembly file
    is provided with `seqin=FILE`. If run with `mapout=FiLE` then an assembly map will be generated from this assembly
    file. If a `*.gaps.tdt` is found (or `gapfile=TDT`), this will be used for gap positions. Otherwise, `rje_seqlist`
    will generate the gaps file and use that. If `telofile=TDT` is provided, telomere positions for the contigs will be
    loaded in. (These must match elements of the assembly map.) If that file is missing, the assembly will be broken into
    contigs and telomere prediction run on the contigs to generate.

    If `mapin=FILE` is given, then GenomeTweak will generate a new assembly (`seqout=FILE` [`$BASEFILE.tweak.fasta]`)
    based on the map file and extracting sequences from the input assembly.

    The map format is as follows:

    ||NewName Description>>SeqName:Start-End:Strand|GapLen| ... |SeqName:Start-End:Strand<<

    Where there is a full-length sequence, Start-End is not required:

    ||NewName Description>>SeqName|~GapLen~| ... |SeqName<<

    Gaps can have zero length.

    If telomere predictions have been loaded from a table (must match SeqName exactly) then 5' and 3' telomeres will be
    annotated in the file with `{` and `}`:

    ||NewName Description>>{SeqName:Start-End:Strand|~GapLen~| ... |SeqName:Start-End:Strand}<<



Commandline:
    ### ~ Main Telociraptor run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly [None]
    basefile=FILE   : Root of output file names [$SEQINBASE]
    tweak=T/F       : Whether to execute GenomeTweak pipeline [True]
    chromsyn=T/F    : Whether to execute ChromSyn preparation (gaps.tdt and telonull telomere prediction) [True]
    dochtml=T/F     : Generate HTML Telociraptor documentation (*.docs.html) instead of main run [False]
    ### ~ Telomere options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    telomeres=T/F   : Whether to generate telomere predictions [True]
    telonull=T/F    : Whether to output sequences without telomeres to telomere table [False]
    telofile=TDT    : Delimited file of `seqname seqlen tel5 tel3` [$SEQINBASE.telomeres.tdt]
    telofwd=X       : Regex for 5' telomere sequence search [C{2,4}T{1,2}A{1,3}]
    telorev=X       : Regex for 5' telomere sequence search [T{1,3}A{1,2}G{2,4}]
    telosize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
    teloperc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
    ### ~ Genome Tweak run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqout=FILE     : Output sequence filename [$BASEFILE.tweak.fasta]
    mapin=FILE      : Text file input for genome assembly map [None]
    mapout=FILE     : Text file output for genome assembly map [$BASEFILE.ctgmap.txt]
    gapfile=TDT     : Delimited file of `seqname start end seqlen gaplen` [$SEQINBASE.gaps.tdt]
    autofix=T/F     : Whether to try to fix terminal inversions based on telomere predictions [True]
    fixout=FILE     : Text file output for auto-fixed genome assembly map [$BASEFILE.tweak.txt]
    invert=LIST     : List of contigs/regions to invert (in order) []
    invlimit=NUM    : Limit inversion distance to within X% (<=100) or Xbp (>100) of chromosome termini [25]
    trimlimit=NUM   : Limit trimming distance to within X% (<=100) or Xbp (>100) of chromosome termini [5]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, re, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_db, rje_readcore, rje_rmd, rje_seqlist, rje_sequence
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial working version with basic map to/from fasta functionality.
    # 0.1.0 - Added autofix=T/F and invert=LIST options.
    # 0.2.0 - Renamed Telociraptor, added standalone telomere prediction, telonull output and chromsyn mode.
    # 0.3.0 - Modified gap formatting to avoid issues with pure number sequence names.
    # 0.4.0 - Upgraded invert=LIST to List of contigs/regions to invert (in order). Added /../ inversion formatting.
    # 0.4.1 - Fixed a bug when sequences have descriptions.
    # 0.5.0 - Added limit for end-trimming and split fixlimit into invlimit and trimlimit.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [N] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Update to generate and recognise SynBad map formats.
    # [ ] : Add ability to read in Gap types from input table and exclude from the map.
    # [Y] : Add option to break into contigs and predict telomeres then read in contig telomere file.
    # [Y] : Add autofix option to invert ends with telomeres.
    # [Y] : Add invert=LIST option to invert individual contigs. (More complex inversions will need manual edits.)
    # [Y] : Check and fix the seqout=FILE setting.
    # [Y] : Rename Telociraptor and add telonull=T/F setting.
    # [ ] : Remove unnecessary duplication of sequence summary prior to and then during gap generation.
    # [Y] : Add a gap character to avoid conflicts with pure number sequence names.
    # [Y] : Extend inversions to multiple sequences X..Y and execute in order given. Recognise /../ inversions.
    # [ ] : Add gap normalisation to correct gap lengths.
    # [ ] : Add summary of fixes at the end of the run, along with change in telomere counts.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('Telociraptor', '0.5.0', 'April 2023', '2023')
    description = 'Telomere Prediction and Genome Assembly Editing Tool'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program when called from commandline.
    '''
    Basic Setup of Program when called from commandline:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:### ~ [1] ~ Initial Command Setup & Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = makeInfo()                                   # Sets up Info object with program details
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2 
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: GenomeTweak Class                                                                                       #
#########################################################################################################################
class GenomeTweak(rje_readcore.ReadCore):
    '''
    GenomeTweak Class. Author: Rich Edwards (2023).

    Str:str
    - FixOut=FILE     : Text file output for auto-fixed genome assembly map [$BASEFILE.tweak.txt]
    - GapFile=TDT     : Delimited file of `seqname start end seqlen gaplen` [$SEQINBASE.gaps.tdt]
    - MapIn=FILE      : Text file input for genome assembly map [None]
    - MapOut=FILE     : Text file output for genome assembly map [None]
    - SeqIn=FILE      : Input sequence assembly [None]
    - SeqOut=FILE     : Output sequence filename [$BASEFILE.tweak.fasta]
    - TeloFile=TDT    : Delimited file of `seqname seqlen tel5 tel3` [$SEQINBASE.telomeres.tdt]
    - TeloFwd=X       : Regex for 5' telomere sequence search [C{2,4}T{1,2}A{1,3}]
    - TeloRev=X       : Regex for 5' telomere sequence search [T{1,3}A{1,2}G{2,4}]

    Bool:boolean
    - AutoFix=T/F     : Whether to try to fix terminal inversions based on telomere predictions [True]
    - ChromSyn=T/F    : Whether to execute ChromSyn preparation (gaps.tdt and telonull telomere prediction) [True]
    - Telomeres=T/F   : Whether to include processing of telomeres [True]
    - TeloNull=T/F    : Whether to output sequences without telomeres to telomere table [False]
    - Tweak=T/F       : Whether to execute GenomeTweak pipeline [True]

    Int:integer
    - Invlimit=NUM    : Limit inversion distance to within X% (<=100) or Xbp (>100) of chromosome termini [25]
    - TeloSize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
    - TrimLimit=NUM   : Limit trimming distance to within X% (<=100) or Xbp (>100) of chromosome termini [5]

    Num:float
    - TeloPerc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]

    File:file handles with matching str filenames
    
    List:list
    - Invert=LIST     : List of contigs to invert []

    Dict:dictionary    

    Obj:RJE_Objects
    - SeqIn = rje_seqlist.SeqList object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['FixOut','GapFile','MapIn','MapOut','SeqIn','SeqOut','TeloFile','TeloFwd','TeloRev']
        self.boollist = ['AutoFix','ChromSyn','DocHTML','Telomeres','TeloNull','Tweak']
        self.intlist = ['InvLimit','TeloSize','TrimLimit']
        self.numlist = ['TeloPerc']
        self.filelist = []
        self.listlist = ['Invert']
        self.dictlist = []
        self.objlist = ['SeqIn']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'TeloFwd':'C{2,4}T{1,2}A{1,3}','TeloRev':'T{1,3}A{1,2}G{2,4}'})
        self.setBool({'AutoFix':True,'ChromSyn':True,'Telomeres':True,'TeloNull':False,'Tweak':True})
        self.setInt({'InvLimit':25,'TeloSize':50,'TrimLimit':5})
        self.setNum({'TeloPerc':50})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['FixOut','GapFile','MapIn','MapOut','SeqIn','SeqOut','TeloFile','TeloFwd','TeloRev'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['AutoFix','ChromSyn','DocHTML','Telomeres','TeloNull','Tweak'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['InvLimit','TeloSize','TrimLimit'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                self._cmdReadList(cmd,'perc',['TeloPerc'])  # Percentage, converts to 1-100 scale.
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Invert'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
#i# seqinObj(self,summarise=False): ### Returns the a SeqList object for the SeqIn file
#i# Also checks for pipe characters in sequence names, which would break things!
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # Telociraptor: Telomere Prediction and Genome Assembly Editing Tool

        Telociraptor is a renamed version of GenomeTweak that also has the capacity for standalone telomere prediction (as
        executed by Diploidocus) along with a combined `chromsyn=T` mode for generating the "telonull" telomere output and
        assembly gaps table used as input by [ChromSyn](https://github.com/slimsuite/chromsyn). The original GenomeTweak
        function is now invoked with `tweak=T`. If running in `tweak` mode, `telomeres` and/or `chromsyn` output will be
        executed on the `seqout=FILE` if generated, which will fixed version of the assembly if `autofix=T`.

        ### ~ Telomere finding [telomeres=T] ~ ###

        Telociraptor performs a regex-based search for Telomeres, based on [FindTelomeres](https://github.com/JanaSperschneider/FindTelomeres).
        By default, this looks for a canonical telomere motif of TTAGGG/CCCTAA, allowing for some variation. Telociraptor
        searches for a forward telomere regex sequence of `C{2,4}T{1,2}A{1,3}` at the 5' end, and a reverse sequence at the
        3' end of `T{1,3}A{1,2}G{2,4}`. These can be set with `telofwd=X` and `telorev=X`. For each sequence, Telociraptor
        trims off any trailing Ns and then searches for telomere-like sequences at sequence ends. For each sequence, the
        presence/absence and length of trimming are reported for the 5' end (tel5 and trim5) and 3' end (tel3 and trim3),
        along with the total percentage telomeric sequence (TelPerc).

        Telomeres are marked if at least 50% (`teloperc=PERC`) of the terminal 50 bp (`telosize=INT`) matches the appropriate
        regex. If either end contains a telomere, the total percentage of the sequence matching either regex is calculated as
        TelPerc. Note that this number neither restricts matches to the termini, not includes sequences within predicted
        telomeres that do not match the regex. By default, only sequences with telomeres are output to the `*.telomeres.tdt`
        output, but switching `telonull=T` will output all sequences. This can be useful if you also need a table of sequence
        lengths, and is the recommended input for [ChromSyn](https://github.com/slimsuite/chromsyn).

        ### ~ GenomeTweak mode [tweak=T] ~ ###

        GenomeTweak is a simple tool designed to help manual curation of genome assembly scaffolding. A genome assembly file
        is provided with `seqin=FILE`. If run with `mapout=FiLE` then an assembly map will be generated from this assembly
        file. If a `*.gaps.tdt` is found (or `gapfile=TDT`), this will be used for gap positions. Otherwise, `rje_seqlist`
        will generate the gaps file and use that. If `telofile=TDT` is provided, telomere positions for the contigs will be
        loaded in. (These must match elements of the assembly map.) If that file is missing, the assembly will be broken into
        contigs and telomere prediction run on the contigs to generate.

        If `mapin=FILE` is given, then GenomeTweak will generate a new assembly (`seqout=FILE` [`$BASEFILE.tweak.fasta]`)
        based on the map file and extracting sequences from the input assembly.

        The map format is as follows:

        ```
        ||NewName Description>>SeqName:Start-End:Strand|~GapLen~| ... |SeqName:Start-End:Strand<<
        ```

        Where there is a full-length sequence, Start-End is not required:

        ```
        ||NewName Description>>SeqName|~GapLen~| ... |SeqName<<
        ```

        Gaps can have zero length.

        If telomere predictions have been loaded from a table (must match SeqName exactly) then 5' and 3' telomeres will be
        annotated in the file with `{` and `}`:

        ```
        ||NewName Description>>{SeqName:Start-End:Strand|~GapLen~| ... |SeqName:Start-End:Strand}<<
        ```

        ## Citation

        Telociraptor has not yet been published. Please cite github in the meantime.

        ---

        # Running Telociraptor

        Telociraptor is written in Python 2.x or Python 3.x and can be run directly from the commandline:

            python $CODEPATH/telociraptor.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [Telociraptor git repo](https://github.com/slimsuite/telociraptor), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [Telociraptor git repo](https://github.com/slimsuite/telociraptor)
        for running on example data.

        ## Dependencies

        The main Telociraptor functions do not currently have any dependencies. To generate documentation with `dochtml`,
        R will need to be installed and a pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        For Telociraptor documentation, run with `dochtml=T` and read the `*.docs.html` file generated.

        ## Commandline options

        A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
        [SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
        use commandline options, including setting default values with **INI files**.

        ```
        ### ~ Main Telociraptor run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FILE      : Input sequence assembly [None]
        basefile=FILE   : Root of output file names [$SEQINBASE]
        tweak=T/F       : Whether to execute GenomeTweak pipeline [False]
        chromsyn=T/F    : Whether to execute ChromSyn preparation (gaps.tdt and telonull telomere prediction) [False]
        dochtml=T/F     : Generate HTML Telociraptor documentation (*.docs.html) instead of main run [False]
        ### ~ Telomere options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        telomeres=T/F   : Whether to generate telomere predictions [True]
        telonull=T/F    : Whether to output sequences without telomeres to telomere table [False]
        telofile=TDT    : Delimited file of `seqname seqlen tel5 tel3` [$SEQINBASE.telomeres.tdt]
        telofwd=X       : Regex for 5' telomere sequence search [C{2,4}T{1,2}A{1,3}]
        telorev=X       : Regex for 5' telomere sequence search [T{1,3}A{1,2}G{2,4}]
        telosize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
        teloperc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
        ### ~ Genome Tweak run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqout=FILE     : Output sequence filename [$BASEFILE.tweak.fasta]
        mapin=FILE      : Text file input for genome assembly map [None]
        mapout=FILE     : Text file output for genome assembly map [$BASEFILE.ctgmap.txt]
        gapfile=TDT     : Delimited file of `seqname start end seqlen gaplen` [$SEQINBASE.gaps.tdt]
        autofix=T/F     : Whether to try to fix terminal inversions based on telomere predictions [True]
        fixout=FILE     : Text file output for auto-fixed genome assembly map [$BASEFILE.tweak.txt]
        invert=LIST     : List of contigs/regions to invert (in order) []
        invlimit=NUM    : Limit inversion distance to within X% (<=100) or Xbp (>100) of chromosome termini [25]
        trimlimit=NUM   : Limit trimming distance to within X% (<=100) or Xbp (>100) of chromosome termini [5]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        # Telociraptor workflow and options

        Details of the Telociraptor workflow will be added with time. Please contact the author or raise an issue on
        GitHub if you have any questions.


        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return rje_rmd.docHTML(self)
            if not self.setup():
                self.printLog('#FAIL','{0} setup failed.'.format(self.prog()))
                return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.seqinObj().setStr({'SeqDictType':'short'})
            ## ~ [2a] GenomeTweak Mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Tweak'):
                self.headLog('TELOCIRAPTOR GENOMETWEAK PROCESSING', line='=')
                if self.getStrLC('MapOut'):
                    ctgmap = self.makeMap()
                    if self.getBool('AutoFix'):
                        self.setStr({'MapIn': self.getStr('MapOut')})
                    else:
                        if not ctgmap: return False
                if self.getStrLC('MapIn'):
                    if self.getBool('AutoFix'):
                        if not self.autoFix():
                            return False
                        self.setStr({'MapIn':self.getStr('FixOut')})
                    if not self.mapToSeq(): return False
            ## ~ [2b] ChromSyn/Telomere mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Tweak') and self.getStrLC('SeqOut'):
                self.setStr({'SeqIn':self.getStr('SeqOut')})
                if self.getBool('ChromSyn') or self.getBool('Telomeres'):
                    self.printLog('#SEQIN','Sequence input updated: {0}'.format(self.getStr('SeqIn')))
            seqbase = rje.baseFile(self.getStr('SeqIn'),strip_path=True)
            if self.getBool('ChromSyn'):
                self.headLog('CHROMSYN GAP GENERATION', line='=')
                self.obj['SeqIn'] = None
                self.setStr({'GapFile':'None'})
                if self.db('gaps'): self.db('gaps').rename('ingaps')
                self.loadGaps()
                #self.seqinObj(summarise=True)
            if self.getBool('ChromSyn') or self.getBool('Telomeres'):
                self.headLog('TELOCIRAPTOR TELOMERE PREDICTION', line='=')
                if self.db('telomeres'): self.db('telomeres').rename('ctgtelomeres')
                return self.findTelomeres(telbase=seqbase,keepnull=self.getBool('TeloNull'))
            elif not self.getBool('Tweak'):
                self.printLog('#RUN','Telociraptor run with telomeres=F tweak=F chromsyn=F. No output.')
                return False
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('SeqIn'): raise IOError('seqin=FILE required.')
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if not self.baseFile(return_none=''):
                if self.getStrLC('SeqIn'): self.baseFile(rje.baseFile(self.getStr('SeqIn'),strip_path=True))
                else: self.baseFile('telociraptor')
            self.printLog('#BASE','Output file basename: %s' % self.baseFile())
            ## ~ [1a] Check/summaries run modes and input files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#MODE','GenomeTweak assembly map/tweaking (tweak=T/F): {0}'.format(self.getBool('Tweak')))
            self.printLog('#MODE','Telomere prediction  (telomeres=T/F): {0}'.format(self.getBool('Telomeres')))
            self.printLog('#MODE','ChromSyn gap and telomere table generation (chromsyn=T/F): {0}'.format(self.getBool('Telomeres')))
            for infile in ['SeqIn','MapIn','GapFile','TeloFile']:
                if self.getStrLC(infile) and not rje.exists(self.getStr(infile)):
                    self.warnLog('{0}={1} set: file not found!'.format(infile,self.getStr(infile)))
            ### ~ [2] Load Sequences and Generate Required Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Tweak'):
                if not self.getStrLC('SeqOut') and self.getStrLC('MapIn'):
                    self.setStr({'SeqOut': '{0}.tweak.fasta'.format(self.baseFile())})
                if not self.getStrLC('MapOut') and not self.getStrLC('MapIn'):
                    self.setStr({'MapOut':'{0}.ctgmap.txt'.format(self.baseFile())})
                if not self.getStrLC('FixOut') and self.getBool('AutoFix'):
                    self.setStr({'FixOut':'{0}.tweak.txt'.format(self.baseFile())})
                    if not self.getStrLC('SeqOut'):
                        self.setStr({'SeqOut': '{0}.tweak.fasta'.format(self.baseFile())})
                if self.getStrLC('MapOut'): ### Gaps needed
                    if self.getStrLC('MapIn'):
                        self.printLog('#MAPIN','Set mapin=None as mapout given.')
                    self.setStr({'MapIn':False})
                    if not self.loadGaps(): return False
                    if not self.loadTelomeres(): return False
            if not self.seqinObj(): return False
            ### ~ [3] Check sequence names for additional issues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin = self.seqinObj()
            for seq in seqin.seqs():
                if rje.matchExp('^~(\d+)~$', seqin.shortName(seq)):
                    self.ValueError('"~$NUMBER~" ({0}) sequence names will cause issues. Please rename and try again.'.format(rje.matchExp('^~(\d+)~$', seqin.shortName(seq))[0]))
            seqin.makeSeqNameDic('short')
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def loadGaps(self): ### Loads assembly gap data, generating if required.
        '''
        Loads assembly gap data, generating if required.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if self.getStrLC('GapFile'):
                gapfile = self.getStr('GapFile')
            else:
                seqbase =  rje.baseFile(self.getStr('SeqIn'),strip_path=True)
                gapfile = '%s.gaps.tdt' % seqbase
            ## ~ [1a] ~ Check or create gapstats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not rje.checkForFiles(filelist=[gapfile],log=self.log):
                self.cmd_list.append('gapstats')    # Try to run automatically if possible
                self.cmd_list.append('seqin={0}'.format(self.getStr('SeqIn')))    # Try to run automatically if possible
                seqin = self.seqinObj(summarise=True)
                if not rje.checkForFiles(filelist=[gapfile],log=self.log):
                    seqin.setBool({'Raw':False,'GapStats':True,'DNA':True})
                    seqin.str['SeqType'] = 'dna'
                    seqin.summarise()
            if not rje.checkForFiles(filelist=[gapfile],log=None):
                raise ValueError('Problem finding/generating {0}'.format(gapfile))

            ### ~ [2] ~ Load Gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdb = db.addTable(gapfile,mainkeys=['seqname','start','end'],name='gaps',ignore=[],expect=True)
            gdb.dataFormat({'seqlen':'int','start':'int','end':'int'})
            return True
        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    ### <3> ### GenomeTweak Methods                                                                                #
#########################################################################################################################
    def makeMap(self): ### Generates assembly map from gap and optionally telomere tables.
        '''
        Generates assembly map from gap and optionally telomere tables.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gapdb = self.db('gaps')      # ['seqname','start','end']
            gapdb.index('seqname')
            teldb = self.db('telomeres') # ['SeqName']
            seqin = self.seqinObj()
            ### ~ [2] Make contig map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            maplist = []
            sx = 0.0
            stot = seqin.seqNum()
            for seq in seqin.seqs():
                self.progLog('\r#CONTIG','Generating sequence contigs for map: %.1f%%' % (sx/stot)); sx += 100.0
                sname = seqin.shortName(seq)
                seqlen = seqin.seqLen(seq)
                seqi = 1
                # Generate Contigs
                seqmap = []
                for entry in gapdb.indexEntries('seqname',sname):
                    seqj = entry['start'] - 1
                    seqmap.append('{0}:{1}-{2}:+'.format(sname,seqi,seqj))
                    seqmap.append('~{0}~'.format(entry['end']-entry['start']+1))
                    seqi = entry['end'] + 1
                if seqi == 1:
                    seqmap.append(sname)
                else:
                    seqmap.append('{0}:{1}-{2}:+'.format(sname, seqi, seqlen))
                # Add telomeres
                for i, contig in enumerate(seqmap):
                    ckey = '.'.join(contig.split(':')[:2])
                    tentry = teldb.data(ckey)
                    if tentry and tentry['Tel5']:
                        seqmap[i] = '{' + '{0}'.format(seqmap[i])
                    if tentry and tentry['Tel3']:
                        seqmap[i] = '{0}'.format(seqmap[i]) + '}'
                # Add to maplist
                maplist.append('||{0}>>{1}<<'.format(seqin.seqName(seq),'|'.join(seqmap)))
            self.printLog('\r#CONTIG', 'Generation of sequence contigs for map complete.')
            ### ~ [3] Output contig map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapout = self.getStr('MapOut')
            rje.backup(self,mapout)
            open(mapout,'w').write('\n'.join(maplist+['']))
            self.printLog('\r#MAPOUT', 'Sequence contig map output to: {0}.'.format(mapout))
            return True
        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    def mapElementData(self,mapel): ### Returns a list of [SeqName, seqi, seqj, Rev, Tel5, Tel3]
        '''
        Returns a list of [SeqName, seqi, seqj, Rev, Tel5, Tel3].
        '''
        seqin = self.seqinObj()
        seqdict = seqin.seqNameDic()
        eldata = ['', 1, 0, None, False, False]
        if mapel[:1] == '{': mapel = mapel[1:]; eldata[4] = True
        if mapel[-1:] == '}': mapel = mapel[:-1]; eldata[5] = True
        if rje.matchExp('^~(\d+)~$',mapel):
            eldata[0] = 'Gap'
            eldata[2] = int(rje.matchExp('^~(\d+)~$', mapel)[0])
        elif ':' in mapel or mapel in seqdict.keys():
            mdata = mapel.split(':')
            if len(mdata) == 3:
                [sname, pos, strand] = mdata
                [seqi, seqj] = pos.split('-')
            else:
                if len(mdata) == 2:
                    [sname, strand] = mdata
                else:
                    sname = mdata[0]
                    strand = '+'
                [seqi, seqj] = [1, seqin.seqLen(seqdict[sname])]
            seqi = int(seqi)
            seqj = int(seqj)
            eldata[0] = sname
            eldata[1] = seqi
            eldata[2] = seqj
            eldata[3] = strand == '-'
        elif mapel.isdigit():
            eldata[0] = 'Gap'
            eldata[2] = int(mapel)
        return eldata
#########################################################################################################################
    def mapLength(self,seqmap): ### Returns the total length of a seqmap
        '''
        Returns the total length of a seqmap
        '''
        maplen = 0
        for mapel in seqmap.split('|'):
            mdata = self.mapElementData(mapel)
            maplen += (mdata[2] - mdata[1] + 1)
        return maplen
#########################################################################################################################
    def invertElement(self,mapel):  ### Inverts a map element
        '''
        Inverts a map element.
        '''
        if rje.matchExp('^~(\d+)~$', mapel): return mapel
        maptel = [False,False]
        if mapel[:1] == '{': mapel = mapel[1:]; maptel[0] = True
        if mapel[-1:] == '}': mapel = mapel[:-1]; maptel[1] = True
        mdata = mapel.split(':')
        newseqmap = ''
        if maptel[1]: newseqmap = newseqmap + '{'
        if mdata[-1] == '+': mdata[-1] = '-'
        elif mdata[-1] == '-': mdata[-1] = '+'
        else: mdata.append('-')
        newseqmap = newseqmap + ':'.join(mdata)
        if maptel[0]: newseqmap = newseqmap + '}'
        return newseqmap
#########################################################################################################################
    def invertRegion(self,seqmap):  ### Reverse order and inverts map elements
        '''
        Reverse order and inverts map elements.
        '''
        invchunk = seqmap.split('|')
        invchunk.reverse()
        newchunk = []
        for i, mapel in enumerate(invchunk):
            newchunk.append(self.invertElement(mapel))
        return '|'.join(newchunk)
#########################################################################################################################
    def autoFix(self): ### Loads the MapIn file and checks it using fixlimit and invert=LIST.
        '''
        Loads the MapIn file and checks it using fixlimit and invert=LIST.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            maplist = open(self.getStr('MapIn'),'r').read()
            maplist = ''.join(maplist.split('\n'))
            seqin = self.seqinObj()
            seqdict = seqin.seqNameDic()
            ## ~ [1a] Split on end of sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            maplist = maplist.split('<<')
            if not maplist[-1].startswith('||'):  maplist = maplist[:-1]
            self.printLog('#CTGMAP','{0} sequences read from {1}'.format(rje.iLen(maplist),self.getStr('MapIn')))

            ### ~ [2] ~ Generate new map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            modified = False
            fixlist = []
            ## ~ [2a] Inversions from /../ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ix = 0
            seqmap = '<<'.join(maplist)
            while rje.matchExp('/([^/]+)/',seqmap):
                invhit = rje.matchExp('/([^/]+)/',seqmap)
                invmap = self.invertRegion(invhit[0])
                self.printLog('#INVERT', '{0} -> {1}'.format(invhit[0], invmap))
                seqmap = seqmap.replace('/{0}/'.format(invhit[0]), invmap)
                modified = True
                ix += 1
            if ix:
                self.printLog('#INVERT','{0} inversions from /../ regions in contig map.'.format(ix))

            ## ~ [2b] Inversions from list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ix = 0
            for invreg in self.list['Invert']:
                invreg = invreg.replace('+','\+')
                if '..' in invreg:
                    invreg = '\S+'.join(invreg.split('..'))
                if rje.matchExp('^(\S+)\.(\d+-\d+)$',invreg):
                    [ctg,ctgreg] = rje.matchExp('^(\S+)\.(\d+-\d+)$',invreg)
                    altreg = ':'.join([ctg,ctgreg])
                    self.printLog('#INVERT', 'Modifying {0} to {1} for inversion'.format(invreg,altreg))
                    invreg = altreg
                invhit = rje.matchExp('[\|>](\{?%s\}?)[\|<]' % invreg,seqmap)
                invhit2 = rje.matchExp('[\|>](\{?%s:[\+-]\}?)[\|<]' % invreg,seqmap)
                if invhit2: invhit = invhit2
                if invhit:
                    invmap = self.invertRegion(invhit[0])
                    self.printLog('#INVERT', '{0} -> {1}'.format(invhit[0], invmap))
                    seqmap = seqmap.replace(invhit[0], invmap)
                    modified = True; ix += 1
                else:
                    self.warnLog('Cannot find "{0}" in seqmap to invert. Check map. May have been edited by previous inversions.'.format(invreg))
            if self.list['Invert']:
                self.printLog('#INVERT','{0} inversions from input list of {1} invert=LIST regions.'.format(ix,len(self.list['Invert'])))
            maplist = seqmap.split('<<')

            ## ~ [2c] Process each sequence and look for telomere-based edits ~~~~~~~~~~~~~~~~~~~~~ ##
            for newseq in maplist:
                [newname,seqmap] = newseq.split('>>')
                if newname[:2] == '||': newname = newname[2:]
                inmap = seqmap
                seqlen = self.mapLength(seqmap)

                # Inversions
                # seqmap = seqmap.split('|')
                # for i, mapel in enumerate(seqmap):
                #     maptel = [False,False]
                #     if mapel[:1] == '{': mapel = mapel[1:]; maptel[0] = True
                #     if mapel[-1:] == '}': mapel = mapel[:-1]; maptel[1] = False
                #     if ':' in mapel or mapel in seqdict.keys():
                #         mdata = mapel.split(':')
                #         if mdata[0] in self.list['Invert'] or ':'.join(mdata[:2]) in self.list['Invert'] or '.'.join(mdata[:2]) in self.list['Invert']:
                #             newseqmap = ''
                #             if maptel[1]: newseqmap = newseqmap + '{'
                #             if mdata[-1] == '+': mdata[-1] = '-'
                #             elif mdata[-1] == '-': mdata[-1] = '+'
                #             else: mdata.append('-')
                #             newseqmap = newseqmap + ':'.join(mdata)
                #             if maptel[0]: newseqmap = newseqmap + '}'
                #             self.printLog('#INVERT','{0} -> {1}'.format(seqmap[i], newseqmap))
                #             seqmap[i] = newseqmap
                #             ix += 1
                # seqmap = ('|').join(seqmap)

                # Look for 5' inversion
                fixlimit = self.getInt('InvLimit')
                tel5 = seqmap.find('{')    # Whether 5' telomere reached
                tel3 = seqmap.find('}')    # Whether 3' telomere reached
                if tel3 > 0 and (tel3 < tel5 or tel5 < 0) and tel3 < (len(seqmap) - 1):
                    # Need to invert up to tel3 and reverse order
                    invchunk = seqmap[:tel3+1]
                    invlen = self.mapLength(invchunk)
                    lentxt = '({0} of {1})'.format(rje_seqlist.dnaLen(invlen),rje_seqlist.dnaLen(seqlen))
                    if fixlimit <= 100 and (100 * float(invlen) / seqlen) > fixlimit: invchunk = ''
                    if fixlimit > 100 and invlen > fixlimit: invchunk = ''
                    if invchunk:
                        keepchunk = seqmap[tel3 + 1:]
                        newchunk = []
                        invchunk = invchunk.split('|')
                        invchunk.reverse()
                        for i, mapel in enumerate(invchunk):
                            newchunk.append(self.invertElement(mapel))
                            # Old gap representation:
                            #if rje.isEven(i): newchunk.append(self.invertElement(mapel))
                            #else: newchunk.append(mapel)
                        newseqmap = '|'.join(newchunk) + keepchunk
                        self.printLog('#INVTEL','{0}... -> {1}...'.format(seqmap[:tel3 + 10], newseqmap[:tel3 + 10]))
                        seqmap = newseqmap
                        modified = True
                    else:
                        self.printLog('#INVTEL','Note: might want to invert 5\' {2} of {0}: {1}'.format(newname.split()[0],seqmap[tel3 + 1:],lentxt))

                # Look for 3' inversion
                tel5 = seqmap.rfind('{')    # Whether 5' telomere reached
                tel3 = seqmap.rfind('}')    # Whether 3' telomere reached
                if tel5 > 0 and (tel5 > tel3):
                    # Need to invert from tel5 and reverse order
                    invchunk = seqmap[tel5:]
                    invlen = self.mapLength(invchunk)
                    lentxt = '({0} of {1})'.format(rje_seqlist.dnaLen(invlen),rje_seqlist.dnaLen(seqlen))
                    if fixlimit <= 100 and (100 * float(invlen) / seqlen) > fixlimit: invchunk = ''
                    if fixlimit > 100 and invlen > fixlimit: invchunk = ''
                    if invchunk:
                        keepchunk = seqmap[:tel5]
                        newchunk = []
                        invchunk = invchunk.split('|')
                        invchunk.reverse()
                        for i, mapel in enumerate(invchunk):
                            newchunk.append(self.invertElement(mapel))
                            # Old gap representation:
                            #if rje.isEven(i): newchunk.append(self.invertElement(mapel))
                            #else: newchunk.append(mapel)
                        newseqmap = keepchunk + '|'.join(newchunk)
                        teli = max(0,tel5-10)
                        self.printLog('#INVTEL','...{0} -> ...{1}'.format(seqmap[teli:], newseqmap[teli:]))
                        seqmap = newseqmap
                        modified = True
                    else:
                        self.printLog('#INVTEL','Note: might want to invert 3\' {2} of {0}: {1}'.format(newname.split()[0],seqmap[tel5:],lentxt))

                # Look for ends to trim off - within fixlimit=INT
                fixlimit = self.getInt('TrimLimit')
                tel5 = seqmap.rfind('{')    # Whether 5' telomere reached
                tel3 = seqmap.rfind('}')    # Whether 3' telomere reached
                if tel5 > 0 and (tel3 > tel5 or tel3 < 0):
                    trimchunk = seqmap[:tel5+1]
                    trimlen = self.mapLength(trimchunk)
                    lentxt = '({0} of {1})'.format(rje_seqlist.dnaLen(trimlen),rje_seqlist.dnaLen(seqlen))
                    if fixlimit <= 100 and (100 * float(trimlen) / seqlen) > fixlimit: trimchunk = ''
                    if fixlimit > 100 and trimlen > fixlimit: trimchunk = ''
                    if trimchunk:
                        self.printLog('#TRIM','Trim 5\' of {0}: {1}'.format(newname.split()[0],trimchunk))
                        seqmap = seqmap[tel5:]
                        for mapel in trimchunk[:tel5].split('|'):
                            if not mapel: continue
                            mdata = self.mapElementData(mapel)
                            if mdata[0] in seqdict:
                                fixlist.append('||{0}.{1}-{2}>>{3}<<'.format(mdata[0],mdata[1],mdata[2],mapel))
                                modified = True
                    else:
                        self.printLog('#TRIM','Note: might want to trim 5\' {2} of {0}: {1}'.format(newname.split()[0],seqmap[:tel5+1],lentxt))
                tel5 = seqmap.find('{')    # Whether 5' telomere reached
                tel3 = seqmap.find('}')    # Whether 3' telomere reached
                endlist = []
                if tel3 < (len(seqmap) - 1) and tel3 > tel5:
                    trimchunk = seqmap[tel3:]
                    trimlen = self.mapLength(trimchunk)
                    lentxt = '({0} of {1})'.format(rje_seqlist.dnaLen(trimlen),rje_seqlist.dnaLen(seqlen))
                    if fixlimit <= 100 and (100 * float(trimlen) / seqlen) > fixlimit: trimchunk = ''
                    if fixlimit > 100 and trimlen > fixlimit: trimchunk = ''
                    if trimchunk:
                        self.printLog('#TRIM','Trim 3\' of {0}: {1}'.format(newname.split()[0],trimchunk))
                        seqmap = seqmap[:tel3+1]
                        for mapel in trimchunk[1:].split('|'):
                            if not mapel: continue
                            mdata = self.mapElementData(mapel)
                            if mdata[0] in seqdict:
                                endlist.append('||{0}.{1}-{2}>>{3}<<'.format(mdata[0], mdata[1], mdata[2], mapel))
                                modified = True
                    else:
                        self.printLog('#TRIM','Note: might want to trim 3\' {2} of {0}: {1}'.format(newname.split()[0],seqmap[tel3:],lentxt))

                fixlist.append('||{0}>>{1}<<'.format(newname,seqmap))
                if endlist: fixlist = fixlist + endlist

            ### ~ [2] ~ Output new map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mapout = self.getStr('FixOut')
            rje.backup(self, mapout)
            open(mapout, 'w').write('\n'.join(fixlist + ['\n']))
            self.printLog('\r#MAPOUT', 'Tweaked sequence contig map output to: {0}.'.format(mapout))
            if not modified: self.printLog('#NOFIX','Tweaked output is identical to input.')

            return True
        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    def mapToSeq(self): ### Loads a contig map and generates sequences from it.
        '''
        Loads a contig map and generates sequences from it.

        The map format is as follows:

        ||NewName Description>>SeqName:Start-End:Strand|GapLen| ... |SeqName:Start-End:Strand<<

        Where there is a full-length sequence, Start-End is not required:

        ||NewName Description>>SeqName|GapLen| ... |SeqName<<

        Gaps can have zero length.

        If telomere predictions have been loaded from a table (must match SeqName exactly) then 5' and 3' telomeres will be
        annotated in the file with `{` and `}`:

        ||NewName Description>>{SeqName:Start-End:Strand|GapLen| ... |SeqName:Start-End:Strand}<<
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin = self.seqinObj()
            seqdict = seqin.seqNameDic()
            maplist = open(self.getStr('MapIn'),'r').read()
            maplist = ''.join(maplist.split('\n'))
            ## ~ [1a] Split on end of sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            maplist = maplist.split('<<')
            if not maplist[-1].startswith('||'):  maplist = maplist[:-1]
            self.printLog('#CTGMAP','{0} sequences read from {1}'.format(rje.iLen(maplist),self.getStr('MapIn')))

            ### ~ [2] ~ Generate and output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.backup(self,self.getStr('SeqOut'))
            SEQOUT = open(self.getStr('SeqOut'),'w')
            for newseq in maplist:
                [newname,seqmap] = newseq.split('>>')
                if newname[:2] == '||': newname = newname[2:]
                SEQOUT.write('>{0}\n'.format(newname))
                newseq = ''
                cx = 1; fx = 0
                for mapel in seqmap.split('|'):
                    if mapel[:1] == '{': mapel = mapel[1:]
                    if mapel[-1:] == '}': mapel = mapel[:-1]
                    if rje.matchExp('^~(\d+)~$', mapel):
                        newseq = newseq + 'N' * int(rje.matchExp('^~(\d+)~$', mapel)[0])
                        cx += 1
                    elif ':' in mapel or mapel in seqdict.keys():
                        mdata = mapel.split(':')
                        if len(mdata) == 3:
                            [sname,pos,strand] = mdata
                            [seqi,seqj] = pos.split('-')
                        else:
                            if len(mdata) == 2:
                                [sname, strand] = mdata
                            else:
                                sname = mdata[0]
                                strand = '+'
                            [seqi, seqj] = [1, seqin.seqLen(seqdict[sname])]
                        seqi = int(seqi)
                        seqj = int(seqj)
                        ctgseq = seqin.seqSequence(seqdict[sname])[seqi-1:seqj]
                        ctglen = seqj - seqi + 1
                        if strand == '-':
                            ctgseq = rje_sequence.reverseComplement(ctgseq,rna=False)
                        if len(ctgseq) != ctglen:
                            raise ValueError('Sequence length mismatch pulling out {0}'.format(mapel))
                        newseq = newseq + ctgseq
                        fx += 1
                    elif int(mapel) > 0:
                        newseq = newseq + 'N' * int(mapel)
                        cx += 1
                SEQOUT.write('{0}\n'.format(newseq))
                self.printLog('#SCAFF','{0} = {1} ({2} fragments; {3} contigs)'.format(newname,rje_seqlist.dnaLen(len(newseq)),fx,cx))
            SEQOUT.close()
            self.printLog('#SEQOUT','{0} sequences saved to {1}'.format(rje.iLen(maplist),self.getStr('SeqOut')))
            return True
        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    ### <4> ### Telomere methods                                                                                        #
#########################################################################################################################
### ~ Telomere options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
#     telomere=X      : Basic telomere sequence for search [TTAGGG]
#     telosize=INT    : Size of terminal regions (bp) to scan for telomeric repeats [50]
#     teloperc=PERC   : Percentage of telomeric region matching telomeric repeat to call as telomere [50]
# parser.add_argument("-w", "--window", type=int, help="This defines the number of first and last nucleotides that will get scanned for telomeric repeats (default: 50).")
# parser.add_argument("-c", "--cutoff", type=float, help='''A telomere is detected if >= c%% of the first (last) nucleotides are telomeric repeats (default: 50%%).''')
#########################################################################################################################
    def loadTelomeres(self): ### Loads assembly telomere data, generating if required.
        '''
        Loads assembly telomere data, generating if required.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('Telomeres'): return False
            seqin = self.seqinObj()
            seqbase = rje.baseFile(self.getStr('SeqIn'), strip_path=True)
            if self.getStrLC('TeloFile'):
                telfile = self.getStr('TeloFile')
            else:
                telfile = '%s.contigs.telomeres.tdt' % seqbase
            ## ~ [1a] ~ Check or create telomeres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not rje.checkForFiles(filelist=[telfile],log=self.log):
                # Contig File
                cfile = '{0}.contigs.fasta'.format(seqbase)
                if not rje.checkForFiles(filelist=[cfile], log=self.log):
                    seqin.saveSeq(seqfile=cfile, reformat='descaffold', append=False, log=True)
            # Telomeres
            seqfile = '{0}.contigs.fasta'.format(seqbase)
            if not rje.exists(seqfile):
                raise IOError('GenomeTweak findTelomeres needs input assembly contigs')
            return self.findTelomeres(seqfile,telbase='{0}.contigs'.format(seqbase),keepnull=self.getBool('TeloNull'))
        except:
            self.errorLog(self.zen())
        return False
#########################################################################################################################
    def findTelomere(self,sequence):    ### Looks for telomeres in nucleotide sequence using telomere regex
        '''
        Looks for telomeres in nucleotide sequence using telomere regex search of sequence ends. Returns a dictionary of
        whether the ends have telomeres and how much was trimmed off as being Ns (5' and 3').
        Based on https://github.com/JanaSperschneider/FindTelomeres.
        >> sequence:str = DNA sequence to search
        << returns dictionary of {'tel5':T/F,'tel3':T/F,'trim5':INT,'trim3':INT,'tel5len':INT,'tel3len':INT}
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tel_forward, tel_reverse = self.getStrUC('TeloFwd'), self.getStrUC('TeloRev')
            sequence = sequence.upper()
            WINDOW = self.getInt('TeloSize')
            REPEAT_CUTOFF = self.getNum('TeloPerc')
            ## ~ [1a] Terminal N-trimming ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            trim5 = 0
            for index, position in enumerate(sequence):
                if position != 'N':
                    trim5 = index
                    break
            start_of_sequence_withoutNs = trim5
            trim3 = 0
            for index, position in enumerate(reversed(sequence)):
                if position != 'N':
                    trim3 = index
                    break
            end_of_sequence_withoutNs = len(sequence) - trim3

            ### ~ [2] Look for Telomeres ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Look for telomeric repeats at the start of the sequence ~~~~~~~~~~~~~~~~~~~~~ ##
            tel5 = 0    # Keep cycling through bigger windows until it breaks down
            while start_of_sequence_withoutNs < end_of_sequence_withoutNs:
                telomeric_repeats = re.findall(tel_forward, sequence[start_of_sequence_withoutNs:start_of_sequence_withoutNs+WINDOW])
                # Calculate the % of nucleotides that are part of telomeric repeats
                percent_telomeric_repeats_start = 100.0*sum([len(repeat) for repeat in telomeric_repeats])/float(WINDOW)
                # If more than half of nucleotides at the start/end are telomeric repeats
                if percent_telomeric_repeats_start >= REPEAT_CUTOFF:
                    tel5 += 1
                    start_of_sequence_withoutNs += WINDOW
                else:
                    break
            telomere_at_start = tel5 > 0
            ## ~ [2b] Look for telomeric repeats at the end of the sequence ~~~~~~~~~~~~~~~~~~~~~ ##
            tel3 = 0
            while end_of_sequence_withoutNs > trim5:
                telomeric_repeats = re.findall(tel_reverse, sequence[(end_of_sequence_withoutNs-WINDOW):end_of_sequence_withoutNs])
                # Calculate the % of nucleotides that are part of telomeric repeats
                percent_telomeric_repeats_end = 100.0*sum([len(repeat) for repeat in telomeric_repeats])/float(WINDOW)
                if percent_telomeric_repeats_end >= REPEAT_CUTOFF:
                    tel3 += 1
                    end_of_sequence_withoutNs -= WINDOW
                else:
                    break
            telomere_at_end = tel3 > 0

            ## ~ [2c] Calculate total percentage telomeres (does not enforce terminal sequences) ~~ ##
            telperc = 0.0
            if telomere_at_start or telomere_at_end:
                telomeric_repeats = re.findall(tel_forward, sequence) + re.findall(tel_reverse, sequence)
                telperc = 100.0 * sum([len(repeat) for repeat in telomeric_repeats]) / float(len(sequence) - sequence.count('N'))

            #!# Update to be more sophisticated and mark end position
            return {'Tel5':telomere_at_start, 'Tel3':telomere_at_end,
                    'Tel5Len':WINDOW*tel5, 'Tel3Len':WINDOW*tel3,
                    'Trim5':trim5, 'Trim3':trim3, 'TelPerc':telperc}
        except:
            self.errorLog('Diploidocus.findTelomere() error'); raise
#########################################################################################################################
    def findTelomeres(self,seqfile=None,save=True,keepnull=False,telbase=None):
        '''
        ## Telomere finding [runmode=telomere]

        Canonical motif is TTAGGG/CCCTAA, but one might see variation, so this searches for regex based on the code at
        https://github.com/JanaSperschneider/FindTelomeres.


        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if not telbase: telbase = db.baseFile()
            seqbase = rje.baseFile(self.getStr('SeqIn'), strip_path=True)
            #i#teldb = self.db().addEmptyTable('telomeres',['SeqName','SeqLen','Tel5','Tel3','Tel5Len','Tel3Len','Trim5','Trim3','TelPerc'],['SeqName'],log=self.debugging())
            telfile = '{}.telomeres.{}'.format(telbase,rje.delimitExt(db.getStr('Delimit')))
            if not self.force() and rje.checkForFiles(filelist=[telfile],basename='',log=self.log):
                teldb = db.addTable(telfile,name='telomeres',mainkeys=['SeqName'])
                teldb.dataFormat({'SeqLen':'int','Tel5':'bool','Tel3':'bool','Trim5':'int','Trim3':'int','TelPerc':'num'})
                return teldb
            forks = self.getInt('Forks')
            ## ~ [1a] ~ Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not seqfile:
                seqfile = self.getStr('SeqIn')
            if not rje.exists(seqfile):
                raise IOError('Telociraptor findTelomeres input assembly not found: {0}'.format(seqfile))
            seqin = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin={0}'.format(seqfile),'autoload=T','seqmode=file','summarise=F','autofilter=F'])
            ## ~ [1b] ~ Results table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            teldb = self.db().addEmptyTable('telomeres',['SeqName','SeqLen','Tel5','Tel3','Tel5Len','Tel3Len','Trim5','Trim3','TelPerc'],['SeqName'],log=self.debugging())
            telomeres = []  # List of sequences with telomeres
            tel5 = tel3 = telboth = 0

            ### ~ [2] ~ Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tel_forward, tel_reverse = self.getStrUC('TeloFwd'), self.getStrUC('TeloRev')
            self.printLog('#TEL','Forward (5\') telomere sequence: {0}'.format(tel_forward))
            self.printLog('#TEL','Reverse (3\') telomere sequence: {0}'.format(tel_reverse))
            sx = 0.0; stot = seqin.seqNum()
            while seqin.nextSeq():
                self.progLog('\r#TELO','Analysing {} sequences for telomeric repeats: {:.2f}%'.format(rje.iStr(stot),sx/stot)); sx += 100.0
                sname = seqin.shortName()
                sequence = seqin.seqSequence()
                tentry = teldb.addEntry(rje.combineDict({'SeqName':sname,'SeqLen':len(sequence)},self.findTelomere(sequence)))
                # Add reporting in verbose mode?
                if tentry['Tel5'] or tentry['Tel3']:
                    telomeres.append(sname)
                    if tentry['Tel5']: tel5 += 1
                    if tentry['Tel3']: tel3 += 1
                    if tentry['Tel5'] and tentry['Tel3']: telboth += 1
                    if self.v() > 0 and tentry['Tel5']:
                        self.printLog('\r#TELO','5\' Telomeric repeat found in {} ({} terminal Ns)'.format(sname,tentry['Trim5']))
                    if self.v() > 0 and tentry['Tel3']:
                        self.printLog('\r#TELO','3\' Telomeric repeat found in {} ({} terminal Ns)'.format(sname,tentry['Trim3']))
            if telomeres: self.printLog('\r#TELO','Telomeric repeats found in {} of {} sequences: {} at 5\', {} at 3\' ({} both)'.format(rje.iLen(telomeres),rje.iStr(stot),tel5,tel3,telboth))
            else: self.printLog('\r#TELO','Telomeric repeats found in {} of {} sequences.'.format(rje.iLen(telomeres),rje.iStr(stot)))

            ### ~ [3] Save Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not keepnull:
                for ekey in list(teldb.dict['Data'].keys()):
                    if ekey not in telomeres: teldb.dict['Data'].pop(ekey)
            if save: teldb.saveToFile(telfile,sfdict={'TelPerc':4})

            return teldb
        except:
            self.errorLog('GenomeTweak.findTelomeres() error')
            return None
#########################################################################################################################
### End of SECTION II: GenomeTweak Class                                                                                #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: GenomeTweak(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
