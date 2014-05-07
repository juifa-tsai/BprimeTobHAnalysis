#!/bin/tcsh
if ( $1 == "" ) then
	echo "./changeName.csh [list] [input] [output]"
	echo "The list structure should be 'fullName;name;"
	exit	
endif
set list_=`cat $1`
foreach list($list_)
	set root=`echo $list | awk -F ";" '{print $1}'`
	set name=`echo $list | awk -F ";" '{print $2}'`
	if ( -e $2/$root.root ) then
		echo "change $root to $name.root..."
		mv $2/$root.root $3/$name.root
	endif
end

