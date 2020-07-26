CurPath=`pwd`
ParPath=$(dirname "$CurPath")

export PATH=$PATH:$ParPath/bin
export PYTHONPATH=$PYTHONPATH:$ParPath
