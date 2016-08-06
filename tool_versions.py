from localutils import *

from subprocess import CalledProcessError

def shell_output_stripped(cmd):
    return check_output_decoded(cmd, shell=True).strip()

# Determine the versions of various programs used
try:
    CUFFLINKS_VERSION = shell_output_stripped('''cufflinks 2>&1 | grep "cufflinks v"''')
except CalledProcessError:
    CUFFLINKS_VERSION = None

try:
    SAMTOOLS_VERSION = shell_output_stripped('''samtools 2>&1 | perl -lane 'print "samtools v$1" if m/Version: (.*)/' ''')
except CalledProcessError:
    SAMTOOLS_VERSION = None

try:
    BOWTIE1_VERSION = shell_output_stripped('''bowtie --version 2>/dev/null | head -n1''')
except CalledProcessError:
    BOWTIE1_VERSION = None

try:
    BOWTIE2_VERSION = shell_output_stripped('''bowtie2 --version 2>/dev/null | head -n1 | perl -lane 'print "bowtie2 $1" if m/(version .*)/' ''')
except CalledProcessError:
    BOWTIE2_VERSION = None

try:
    TOPHAT2_VERSION = shell_output_stripped('''tophat2 --version 2>/dev/null''')
except CalledProcessError:
    TOPHAT2_VERSION = None

try:
    HISAT2_VERSION = shell_output_stripped('''hisat2 --version 2>/dev/null | head -n1 | perl -lane 'print "hisat2 $1" if m/(version .*)/' ''')
except CalledProcessError:
    HISAT2_VERSION = None

try:
    BWA_VERSION = shell_output_stripped('''bwa 2>&1 | perl -lane 'print "bwa v$1" if m/Version: (.*)/' ''')
except CalledProcessError:
    BWA_VERSION = None

try:
    BBMAP_VERSION = shell_output_stripped('''{BBMAP} --version 2>&1 | grep 'BBMap version' '''.format(BBMAP=BBMAP))
except CalledProcessError:
    BBMAP_VERSION = None

try:
    STAR_VERSION = shell_output_stripped('''STAR --version 2>/dev/null''')
except CalledProcessError:
    STAR_VERSION = None

try:
    SALMON_VERSION = 'salmon ' + shell_output_stripped('''salmon --version 2>&1''')
except CalledProcessError:
    SALMON_VERSION = None

try:
    from rpy2.robjects import r
    from rpy2.rinterface import RRuntimeError
    R_VERSION = ''.join(r('R.version$version.string'))
except RRuntimeError:
    R_VERSION = None

try:
    from rpy2.robjects import r
    from rpy2.rinterface import RRuntimeError
    BIOC_VERSION = 'Bioconductor ' + "".join(r('as.character(BiocInstaller::biocVersion())'))
except RRuntimeError:
    BIOC_VERSION = None
