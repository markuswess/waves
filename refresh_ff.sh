WID=`xdotool search --name "Mozilla Firefox" | head -1`
xdotool windowactivate $WID
xdotool key F5
