from localutils import *

from rpy2 import robjects
from rpy2.rinterface import RRuntimeError

def shell_output_stripped(cmd):
    return check_output_decoded(cmd, shell=True).strip()

# Determine the versions of various programs used
CUFFLINKS_VERSION = shell_output_stripped('''cufflinks 2>&1 | grep "cufflinks v"''')
SAMTOOLS_VERSION = shell_output_stripped('''samtools 2>&1 | perl -lane 'print "samtools v$1" if m/Version: (.*)/' ''')
BOWTIE1_VERSION = shell_output_stripped('''bowtie --version 2>/dev/null | head -n1''')
BOWTIE2_VERSION = shell_output_stripped('''bowtie2 --version 2>/dev/null | head -n1 | perl -lane 'print "bowtie2 $1" if m/(version .*)/' ''')
TOPHAT2_VERSION = shell_output_stripped('''tophat2 --version 2>/dev/null''')
HISAT2_VERSION = shell_output_stripped('''hisat2 --version 2>/dev/null | head -n1 | perl -lane 'print "hisat2 $1" if m/(version .*)/' ''')
BWA_VERSION = shell_output_stripped('''bwa 2>&1 | perl -lane 'print "bwa v$1" if m/Version: (.*)/' ''')
BBMAP_VERSION = shell_output_stripped('''{BBMAP} --version 2>&1 | grep 'BBMap version' '''.format(BBMAP=BBMAP))
STAR_VERSION = shell_output_stripped('''STAR --version 2>/dev/null''')
SALMON_VERSION = 'salmon ' + shell_output_stripped('''salmon --version 2>&1''')

R_VERSION = ''.join(robjects.r('R.version$version.string'))
try:
    BIOC_VERSION = 'Bioconductor ' + "".join(robjects.r('as.character(BiocInstaller::biocVersion())'))
except RRuntimeError:
    BIOC_VERSION = None
