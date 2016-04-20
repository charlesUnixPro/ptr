all : doit

doit: ptr frames process
	rm process/*.pam
	./ptr
	rm out.mpg
	ffmpeg -framerate 30 -i process/6_straight%05d.pam -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4

ptr: ptr.c frames
	clang ptr.c -o ptr -lnetpbm  -O3 -g  -lm

frames:
	mkdir frames
	./1_convert2pam

process:
	mkdir process
