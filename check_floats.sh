
#$Commands defines everything ran on the seadata.fmi.fi.
#$Cmnd is formatted version to ensure seadata's broken bash actually gets it all.

Commands="
echo ------------------------;
echo FILES/DATES:;
echo ------------------------;
last_log=\$(ls *.log -rth |tail -n 1);
last_msg=\$(ls *.msg -rth |tail -n 1);
ls *.log -lrth |tail -n 1;
ls *.msg -lrth |tail -n 1;
echo Last Transferred: ;
cat LastFileTransferred.dat;
echo;
echo Last Processed   :;
cat LastFileProcessed.dat;
echo;
echo ------------------------;
echo PRESSURES:;
echo ------------------------;
grep ParkTerminate \$last_log;
grep ProfileInit \$last_log;
grep ParkPts \$last_msg;
echo ------------------------;
echo LOCATION: ;
echo ------------------------;
grep Fix: \$last_msg
"


Cmnd=$(echo "echo '${Commands}'" | tr -d "\n" )

if [ -z $1 ]  ; 
then
	Float="f7126"
else    
	Float=$1
fi
ssh seadata.fmi.fi -l $Float "${Cmnd} | bash"
#Float="f7087"
#ssh seadata.fmi.fi -l $Float "${Cmnd} | bash"