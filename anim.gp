do for[a=1:1217]{
	pl 'sigmoid_3d.txt' u 1:2:3 every :::a::a w p pt 7 palette
	pause 0.1
}
pause -1
