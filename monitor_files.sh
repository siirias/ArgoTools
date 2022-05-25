old=$(ssh seadata.fmi.fi -l f7126 'ls -lrth --color=none|tail -n 1') ; 
while [ 1 -gt 0 ]; do 
	new=$(ssh seadata.fmi.fi -l f7126 'ls -lrth --color=none|tail -n 1') ; 
	if [ "$new" != "$old" ] ; then
		zenity --notification --text="new file $new" ;
		fi;
	old=$new ; 
	sleep 10;
done