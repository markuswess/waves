WID=`xdotool search --name "Numerical methods for wave-type equations" | head -1`
WIDVI=`xdotool search --name "jupyter-book" | head -1`
xdotool windowactivate $WID
xdotool key F5
xdotool windowactivate $WIDVI
