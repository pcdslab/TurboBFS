TEST="1"
DATA="/lclhome/oarti001/documents/cuBFS/graphData" #update this path

for file in mycielskian4
do
    if [ "$TEST" = "1" ] ; then
	echo $file
	echo TurboBFS_CSC_TDBU nvprof --profile-api-trace none  ./exec $DATA/$file.mtx --w 0 --ug  1  --format  1 --p 1 --r 0 --repet 1 --seq 1 --bbreak 1
	nvprof --profile-api-trace none  ./exec $DATA/$file.mtx  0 1 1 1 0 1 1 1
    fi
done

for file in Tina_AskCal
do
    if [ "$TEST" = "1" ] ; then
	echo $file
	echo TurboBFS_CSC_TDBU nvprof --profile-api-trace none  ./exec $DATA/$file.mtx --w 0 --ug  0 --format  1 --p 1 --r 10 --repet 1 --seq 1 --bbreak 1
	nvprof --profile-api-trace none  ./exec $DATA/$file.mtx  0 0 1 1 10 1 1 1
    fi
done

for file in Ragusa18 cage4
do
    if [ "$TEST" = "1" ] ; then
	echo $file
	echo TurboBFS_CSC_TDBU nvprof --profile-api-trace none  ./exec $DATA/$file.mtx --w 1 --ug  0 --format  1 --p 1 --r 0 --repet 10 --seq 1 --bbreak 1
	nvprof --profile-api-trace none  ./exec $DATA/$file.mtx  1 0 1 1 0 10 1 1
    fi
done

for file in smallworld luxembourg_osm 
do
    if [ "$TEST" = "2" ] ; then
	echo $file
	echo TurboBFS_CSC_TDBU nvprof --profile-api-trace none  ./exec $DATA/$file.mtx --w 0 --ug 1 --format  0 --p 0 --r 0 --repet 100 --seq 1 --bbreak 1
	nvprof --profile-api-trace none  ./exec $DATA/$file.mtx  0 1 0 0 0 100  1 1
    fi
done

for file in  amazon0302 
do
    if [ "$TEST" = "2" ] ; then
	echo $file
	echo TurboBFS_CSC_TDBU nvprof --profile-api-trace none  ./exec $DATA/$file.mtx --w 0 --ug 0 --format  0 --p 0 --r 10 --repet 100 --seq 1 --bbreak 1
	nvprof --profile-api-trace none  ./exec $DATA/$file.mtx  0 0 0 0 10 100  1 1
    fi
done
