#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <jpeglib.h>	// apt install libjpeg-dev
#include <fftw3.h>	// from fftw.org
#include <libconfig.h>	// apt install libconfig-dev

typedef unsigned char BYTE;
typedef unsigned int bool;

#define	CONFIG_FILE	"config"
config_t Config;

// Read from config file
unsigned int	Width_out;	// width of cropped output image
unsigned int	Height_out;	// height of cropped output image
unsigned int	Left_out;	// left column of cropped output image
unsigned int	Top_out;	// top row of cropped output image
unsigned int	Fft_width;	// number of columns of fft window (must be power of 2, >=2)
unsigned int	Fft_height;	// number of rows of fft window (must be power of 2, >=2)
unsigned int	Fft_left;	// left column of FFT window
unsigned int	Fft_top;	// top row of FFT window
unsigned int	Border_size;	// size of gray pixel border
unsigned int	Quality;	// outgoing compression quality (percentage)
int		Verbose;
int		Debug;
int		Plot_correl;	// plot correlation surface
const char	*Align_dir;	// subdirectory to write aligned frames
const char	*Correl_dir;	// subdirectory to write correlation surface image
const char	*Reference_name;// image to use as reference to all the others
double		Red_weight;
double		Blue_weight;
double		Green_weight;
double		Trust;

// Taken from reference frame. All input frames must have identical parameters
unsigned int	Width_in;
unsigned int	Height_in;
unsigned int	Components = 3;	// RGB only, cannot change
unsigned int	Colorspace;	// usually RGB. might be Grayscale?
double		Gamma;

// Derived from reference frame and/or config settings
unsigned int	Row_length;	// computed length of each row (Width_in*Components)
unsigned int	Frame_size;	// computed total size of each frame (Width_in*Height_in*Components)
unsigned int	Fft_size;	// Fft_width * Fft_height
double		All_weight;	// Red_weight + Blue_weight + Green_weight

// stats
unsigned int	Total_frames = 0;
unsigned int	Low_trust = 0;
unsigned int	Big_dx = 0;
unsigned int	Big_dy = 0;

#define	RED_IDX		0
#define	GRN_IDX		1
#define	BLU_IDX		2

// basic colors
#define	COLOR_RED	0xFF0000
#define	COLOR_GREEN	0x00FF00
#define	COLOR_BLUE	0x0000FF
#define	COLOR_YELLOW	0xFFFF00
#define	COLOR_CYAN	0xFF00FF
#define	COLOR_MAGENTA	0x00FFFF
#define	COLOR_GRAY	0x808080
#define	COLOR_WHITE	0xFFFFFF
#define	COLOR_BLACK	0x000000


// FFT
fftw_plan Plan_fwd;			// Forward FFT plan
fftw_plan Plan_inv;			// Inverse FFT plan
#define	RE	0
#define	IM	1
fftw_complex *Ref_fft;			// complex FFT of reference frame
fftw_complex *Correl;			// complex correlation surface
double *Surface;			// real correlation surface

//	print error message and exit
static inline void
fatal(char *format, ...)
{
	va_list args;

	printf("ERROR: ");
	va_start(args, format);
	vprintf(format, args);
	va_end(args);
	printf("\n");
	exit(1);
}

// fopen or fail
static inline FILE *
must_fopen(const char *f,const char *m)
{
	FILE *fp = fopen(f,m);

	if( fp==NULL )
		fatal("Can't fopen %s for mode %s",f,m);
	return fp;
}


// malloc or fail
static inline void *
must_malloc(const size_t size)
{
	void *vp = malloc(size);

	if(vp==NULL)
		fatal("Can't malloc %d bytes",size);
	return vp;
}

// make sure image parameters are all identical
static inline void
image_verify(const unsigned int width, const unsigned int height, const unsigned int components, const unsigned int colorspace, const double gamma)
{
	if( Width_in==0 && Height_in==0 ){	// first frame sets it for the rest
		Width_in	= width;
		Height_in	= height;
		//Components	= components;
		Colorspace	= colorspace;
		Gamma		= gamma;
		Row_length	= Width_in*Components;
		Frame_size	= Height_in*Row_length;
		}
	if( width != Width_in )		fatal("width %u != %u",Width_in,width);
	if( height != Height_in )	fatal("height %u != %u",Height_in,height);
	if( components != Components )	fatal("components %u != %u",Components,components);
	if( colorspace != Colorspace )	fatal("colorspace %u != %u",Colorspace,colorspace);
	if( gamma != Gamma )		fatal("gamma %u != %u",Gamma,gamma);
}

// load image, convert from JPG to uncompressed RGB, return image as a giant string of bytes
BYTE *
load_image (const char *name)
{
	FILE *fp = must_fopen (name, "rb");
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW row_pointer[1];
	BYTE *p;

	cinfo.err = jpeg_std_error (&jerr);
	jpeg_create_decompress (&cinfo);
	jpeg_stdio_src (&cinfo, fp);
	jpeg_read_header (&cinfo, TRUE);
	jpeg_start_decompress (&cinfo);

	image_verify(cinfo.output_width,cinfo.output_height,cinfo.output_components,cinfo.out_color_space,cinfo.output_gamma);

	p = (BYTE *) must_malloc (Frame_size);

	while (cinfo.output_scanline < Height_in) {
		row_pointer[0] = &p[cinfo.output_scanline * Row_length];
		jpeg_read_scanlines (&cinfo, row_pointer, 1);
	}

	jpeg_finish_decompress (&cinfo);
	jpeg_destroy_decompress (&cinfo);
	fclose (fp);
	return p;
}

// crop image according to Top_out, Left_out, Width_out, Height_out
static inline void
crop_image (BYTE *p)
{
	unsigned int row, col;
	unsigned int newstride = Width_out * Components;

	for (row = 0; row < Height_out; row++)
		for (col = 0; col < newstride; col++)
			p[row * newstride + col] = p[((row + Top_out) * Row_length) + ((Left_out * Components) + col)];
}

// write image to file using Width_out, Height_out
void
store_image (char *name, BYTE *p, const unsigned int width, const unsigned int height, const unsigned int components, const unsigned int colorspace, const double gamma)
{
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	FILE *fp;
	JSAMPROW row_pointer[1];

	if (p == NULL)fatal("store_image has no data");

	fp = must_fopen (name, "wb");

	cinfo.err = jpeg_std_error (&jerr);
	jpeg_create_compress (&cinfo);
	jpeg_stdio_dest (&cinfo, fp);

	cinfo.image_width	= width;
	cinfo.image_height	= height;
	cinfo.input_components	= components;
	cinfo.in_color_space	= colorspace;
	cinfo.input_gamma	= gamma;
	jpeg_set_defaults (&cinfo);
	jpeg_set_quality (&cinfo, Quality, TRUE);
	jpeg_start_compress (&cinfo, TRUE);

	while (cinfo.next_scanline < height) {
		row_pointer[0] = p;
		p += width*components;
		jpeg_write_scanlines (&cinfo, row_pointer, 1);
	}

	jpeg_finish_compress (&cinfo);
	jpeg_destroy_compress (&cinfo);
	fclose (fp);
}

// rotate a chunk of memory
static inline void
rotate_mem(BYTE *buf, const unsigned int mid, const unsigned int len)
{
	BYTE *tmp;

	if( mid > len )
		fatal("rotate_mem mid:%u len:%u\n",mid,len);
	if( mid==0 || mid==len )
		return;
	tmp = must_malloc(mid);
	memmove(tmp,buf,mid);		// save the front chunk of buf up to mid
	memmove(buf,buf+mid,len-mid);	// copy back chunk to front
	memmove(buf+(len-mid),tmp,mid);	// put saved chunk at end
	free(tmp);
}

// rotate horizontal, with wraparound
static inline void
rotate_image_horizontal(BYTE *p, int ncols, const unsigned int width, const unsigned int height)
{
	unsigned int row;
	unsigned int row_length = width*Components;

	if(ncols==0)return;
	while(ncols < 0 )
		ncols += Width_in;

	for(row=0; row<height; row++){
		rotate_mem(p,Components*(ncols%width),row_length);
		p += row_length;
		}
}

// rotate vertical, with wraparound
static inline void
rotate_image_vertical(BYTE *p, int nrows, const unsigned int width, const unsigned int height)
{
	unsigned int row_length = width*Components;
	unsigned int frame_size = row_length*height;

	if(nrows==0)return;
	while(nrows < 0 )
		nrows += height;

	rotate_mem(p,row_length*(nrows%height),frame_size);
}

static inline void
add_vertical_line(BYTE *buf, const int col)
{
	unsigned int i;

	buf += Components*(col%Width_in);
	for(i=0;i<Height_in;i++){
		memset(buf,0x80,Components);
		buf += Row_length;
		}
}

static inline void
add_horizontal_line(BYTE *buf, const int row)
{
	buf += Row_length*(row%Height_in);
	memset(buf,0x80,Row_length);
}

static inline void
add_border(BYTE *buf, const unsigned int width)
{
	unsigned int i;

	for(i=0;i<width;i++){
		add_vertical_line(buf,i);
		add_vertical_line(buf,Width_in-i-1);
		add_horizontal_line(buf,i);
		add_horizontal_line(buf,Height_in-i-1);
		}
}

static inline void
draw_dot(BYTE *data, const unsigned int color, const unsigned int x, const unsigned int y)
{
	data += (y%Height_in)*Row_length;
	data += (x%Width_in)*Components;

	data[RED_IDX] = (color>>16)&0xFF;
	data[GRN_IDX] = (color>> 8)&0xFF;
	data[BLU_IDX] = (color>> 0)&0xFF;
}

static inline void
draw_box(BYTE *data, const unsigned int color, const unsigned int top, const unsigned left, const unsigned int width, const unsigned int height)
{
	unsigned int x,y;

	for(x=left; x<(left+width); x++){
		draw_dot(data,color,x,top);
		draw_dot(data,color,x,top+height-1);
	}
	for(y=top; y<(top+height); y++){
		draw_dot(data,color,left,y);
		draw_dot(data,color,left+width-1,y);
	}
}

// convert RGB pixel to grayscale, weighted and scaled [0..1)
static inline double
rgb2real(const BYTE *pixel, const unsigned int row, const unsigned int col)
{
	pixel += (row*Row_length)+(col*Components);
	return ((pixel[RED_IDX]*Red_weight)+(pixel[GRN_IDX]*Green_weight)+(pixel[BLU_IDX]*Blue_weight))/All_weight;
}

// extract FFT window (Fft_top,Fft_left,Fft_width,Fft_height) from frame, convert RGB pixel to real value
static inline void
frame_data2d (const BYTE *data, double *realdata)
{
	unsigned int row,col;

	for (row=0; row<Fft_height; row++)
	for (col=0; col<Fft_width;  col++)
		*realdata++ = rgb2real(data,Fft_top+row,Fft_left+col);
}

// rotate image, add helpful info in debug mode
static inline void
rotate_image(BYTE *data, const int dx, const int dy)
{
	if(Debug){
		add_border(data,Border_size);
		// before alignment
		draw_box(data,COLOR_RED,Fft_top,Fft_left,Fft_width,Fft_height);
		draw_box(data,COLOR_RED,Fft_top-1,Fft_left-1,Fft_width+2,Fft_height+2);
		}
	rotate_image_vertical(data, dy, Width_in, Height_in);			// rotate the image
	rotate_image_horizontal(data, dx, Width_in, Height_in);
	if(Debug){
		// after alignment
		draw_box(data,COLOR_BLUE,Fft_top,Fft_left,Fft_width,Fft_height);
		draw_box(data,COLOR_BLUE,Fft_top-1,Fft_left-1,Fft_width+2,Fft_height+2);
		}
}

// make fft of frame
static inline fftw_complex *
get_plane_fft (const BYTE *data)
{
	fftw_complex *fft = (fftw_complex *)fftw_malloc (sizeof (fftw_complex) * Fft_size);
	double *realdata = (double *)fftw_malloc( sizeof(double) * Fft_size);

	frame_data2d (data, realdata);
	fftw_execute_dft_r2c (Plan_fwd, realdata, fft); // make forward fft of data
	free(realdata);
	return fft;
}

// multiply complex conj. (mult = a * b) (hermetian)
static inline void
mult_conj_data2d (const fftw_complex *a, const fftw_complex *b, fftw_complex *mult)
{
	unsigned int i;

	for (i = 0; i < Fft_size; i++ ) {
		mult[i][RE] = a[i][RE] * b[i][RE] + a[i][IM] * b[i][IM];
		mult[i][IM] = a[i][RE] * b[i][IM] - a[i][IM] * b[i][RE];
		}
}

// convert unsigned value to signed, modulo the 'mod' size
static inline int
signed_modulo(const unsigned int val, const unsigned int mod)
{
	return val<(mod/2) ? val : val-mod;
}

// find global max on real part of correlation surface
// return an estimate of how good the result is (best/average)
static inline double
find_max_correl (const double *correl, int *dx, int *dy)
{
	unsigned int i;
	unsigned int best_idx = 0;
	double best = correl[0];
	double average = correl[0];

	// find max in correlation array
	for( i=1; i<Fft_size; i++){
		average += correl[i];
		if( best < correl[i] ){
			best = correl[i];
			best_idx = i;
			}
		}
	average /= Fft_size;

	*dx = signed_modulo(best_idx%Fft_width,Fft_width);
	*dy = signed_modulo(best_idx/Fft_width,Fft_height);

	return average ? best/average : 0.0;
}

// convert correlation value to gray scale pixel
static inline BYTE
correl_to_gray(const double val, const double worst, const double best)
{
	double step = (best-worst)/256;

	if( val<=worst )return 0x00;
	if( val>=best  )return 0xFF;
	return (BYTE)((val-worst)/step);
}

// convert correlation array to pixel array
// green = good, red = bad
static inline BYTE *
correl_to_pixels(const double *correl)
{
	BYTE *p = must_malloc(Fft_size*Components);
	BYTE *q = p;
	BYTE level;
	unsigned int i;
	double best = correl[0];
	double worst = correl[0];

	// find min and max
	for( i=1; i<Fft_size; i++){
		if( best < correl[i] )
			best = correl[i];
		if( worst > correl[i] )
			worst = correl[i];
		}

	for(i=0;i<Fft_size; i++){
		level = correl_to_gray(correl[i],worst,best);
		q[RED_IDX] = 0xFF-level;
		q[GRN_IDX] = level;
		q[BLU_IDX] = 0;
		q += Components;
		}
	return p;
}

// alignm one frame against the reference
static inline void
align(const char *name)
{
	BYTE *data = load_image(name);			// get the pixels
	fftw_complex *fft = get_plane_fft(data);	// create 2D FFT
	int dx,dy;					// best guess of delta to reference frame
	double trust;					// how much do we trust the answer
	BYTE *p;
	char store_name[1024];
	unsigned int low_trust, big_dx, big_dy;		// stats for this frame

	mult_conj_data2d (fft, Ref_fft, Correl);	// create correlation data (Correl = ref * curr)
	fftw_execute_dft_c2r (Plan_inv, Correl, Surface);	// make inverse fft of correlation data
	trust = find_max_correl (Surface, &dx, &dy);	// dx,dy global motion as max on correlation surface
	rotate_image(data,-dx,-dy);
	crop_image(data);
	sprintf(store_name,"%s/%s",Align_dir,name);
	store_image(store_name,data,Width_out,Height_out,Components,Colorspace,Gamma);
	if( Plot_correl ){
		sprintf(store_name,"%s/%s",Correl_dir,name);
		p = correl_to_pixels(Surface);
		rotate_image_vertical(p,Fft_height/2,Fft_width,Fft_height);
		rotate_image_horizontal(p,Fft_width/2,Fft_width,Fft_height);
		store_image(store_name,p,Fft_width,Fft_height,Components,Colorspace,Gamma);
		free(p);
		}
	low_trust = (trust<Trust) ? 1 : 0;
	big_dx = (abs(dx)>Fft_width/4) ? 1 : 0;
	big_dy = (abs(dy)>Fft_height/4) ? 1 : 0;
	Low_trust += low_trust;
	Big_dx += big_dx;
	Big_dy += big_dy;
	Total_frames++;

	if(Verbose)
		printf("%s %6d %6d %7.3lf%s%s%s\n", name,dx,dy,trust,
			low_trust ? " Trust?":"",
			big_dx    ? " Dx?"   :"",
			big_dy    ? " Dy?"   :"");
	free(fft);
	free(data);
}

// read an unsigned int from the config file, validate against min/max values
void
configure_unsigned_int(char *name, unsigned int *val, const unsigned int minval, const unsigned int maxval)
{
	int v;

	if( config_lookup_int(&Config,name,&v) != CONFIG_TRUE )
		fatal("unknown config name, %s",name);
	if( v<0 )
		fatal("expected non-negative config for %s, got %d",name,v);
	*val = (unsigned int)v;
	if( *val < minval || *val > maxval )
		fatal("invalid range for %s: %u, min:%u max:%u",name,*val,minval,maxval);
}

// read a float from the config file, validate against min/max values
void
configure_double(char *name, double *val, const double minval, const double maxval)
{
	if( config_lookup_float(&Config,name,val) != CONFIG_TRUE )
		fatal("unknown config name, %s",name);
	if( *val < minval || *val > maxval )
		fatal("invalid range for %s: %6.2lf, min:%6.2lf max:%6.2lf",name,*val,minval,maxval);
}

void
configure_string(char *name, const char **val)
{
	const char *p;

	if( config_lookup_string(&Config,name,&p) != CONFIG_TRUE )
		fatal("unknown config name, %s",name);
	*val = must_malloc(strlen(p)+1);	// make copy since config_destroy will release it
	strcpy((char *)*val,p);
}

void
configure()
{
	config_init(&Config);
	if( config_read_file(&Config,CONFIG_FILE) != CONFIG_TRUE )
		fatal("No config file? (%s)",CONFIG_FILE);
	configure_string("reference",&Reference_name);
	configure_unsigned_int("width",&Width_out,1,8192);
	configure_unsigned_int("height",&Height_out,1,8192);
	configure_unsigned_int("left",&Left_out,0,8192);
	configure_unsigned_int("top",&Top_out,0,8192);
	configure_unsigned_int("fft_width",&Fft_width,2,13);
	configure_unsigned_int("fft_height",&Fft_height,2,13);
	configure_unsigned_int("fft_left",&Fft_left,0,8192);
	configure_unsigned_int("fft_top",&Fft_top,0,8192);
	configure_unsigned_int("quality",&Quality,0,100);
	configure_unsigned_int("border",&Border_size,0,4096);
	configure_double("red_weight",&Red_weight,0,1000);
	configure_double("green_weight",&Green_weight,0,1000);
	configure_double("blue_weight",&Blue_weight,0,1000);
	configure_double("trust",&Trust,0,100000);
	configure_string("align_dir",&Align_dir);
	configure_string("correl_dir",&Correl_dir);
	config_lookup_bool(&Config,"verbose",&Verbose);
	config_lookup_bool(&Config,"debug",&Debug);
	config_lookup_bool(&Config,"correlation",&Plot_correl);
	config_destroy(&Config);

	All_weight  = Red_weight + Green_weight + Blue_weight;
}

static inline void
show_settings()
{
	printf("Image parameters\n");
	printf("\tWidth\t\t%u\n",Width_in);
	printf("\tHeight\t\t%u\n",Height_in);
	printf("\tColorspace\t%u\n",Colorspace);
	printf("\tGamma\t\t%5.2lf\n",Gamma);
	printf("\tRow_length\t%u\n",Row_length);
	printf("\tFrame_size\t%u\n",Frame_size);
	printf("Config parameters\n");
	printf("\tFft_width\t%u\n",Fft_width);
	printf("\tFft_height\t%u\n",Fft_height);
	printf("\tFft_left\t%u\n",Fft_left);
	printf("\tFft_top\t\t%u\n",Fft_top);
	printf("\tWidth_out\t%u\n",Width_out);
	printf("\tHeight_out\t%u\n",Height_out);
	printf("\tLeft_out\t%u\n",Left_out);
	printf("\tTop_out\t\t%u\n",Top_out);
	printf("\tBorder_size\t%u\n",Border_size);
	printf("\tRed_weight\t\t%10.3lf\n",Red_weight);
	printf("\tGreen_weight\t\t%10.3lf\n",Green_weight);
	printf("\tBlue_weight\t\t%10.3lf\n",Blue_weight);
	printf("\tTrust\t\t%5.2lf\n",Trust);
	printf("\tQuality\t\t%u\n",Quality);
	printf("\tVerbose\t\t%s\n",Verbose?"true":"false");
	printf("\tDebug\t\t%s\n",Debug?"true":"false");
	printf("\tCorrelation\t%s\n",Plot_correl?"true":"false");
	printf("\tAlign_dir\t%s\n",Align_dir);
	printf("\tCorrel_dir\t%s\n",Correl_dir);
	printf("\tReference\t%s\n",Reference_name);
}

static inline void
fft_setup()
{
	// adjust FFT from log2 to actual
	Fft_width	= 1<<Fft_width;
	Fft_height	= 1<<Fft_height;

	Fft_size	= Fft_width * Fft_height;	// derived value

	// memory for correlation matrices
	Correl = (fftw_complex *) fftw_malloc (sizeof (fftw_complex) * Fft_size );
	Surface = (double *)fftw_malloc( sizeof(double) * Fft_size );

	// create FFTW plans
	Plan_fwd = fftw_plan_dft_r2c_2d (Fft_height, Fft_width, Surface, Correl, FFTW_ESTIMATE);
	Plan_inv = fftw_plan_dft_c2r_2d (Fft_height, Fft_width, Correl, Surface, FFTW_ESTIMATE);
}

// read the reference frame and compute its FFT
static inline void
reference_setup()
{
	BYTE *data = load_image(Reference_name);	// establish incoming frame parameters (all must be identical)
	Ref_fft = get_plane_fft(data);
	free(data);
}

// finalize globals and look for problems
static inline void
check_config()
{
	if(Verbose)
		show_settings();

	// individual range checks have passed.  check combinations for validity
	if( (Fft_width+Fft_left) > Width_in )fatal("FFT width(%d)+left(%d) > %d",Fft_width,Fft_left,Width_in);
	if( (Fft_height+Fft_top) > Height_in )fatal("FFT height(%d)+top(%d) > %d",Fft_height,Fft_top,Height_in);
	if( (Width_out+Left_out) > Width_in )fatal("out width(%d)+left(%d) > %d",Width_out,Left_out,Width_in);
	if( (Height_out+Top_out) > Height_in )fatal("out height(%d)+top(%d) > %d",Height_out,Top_out,Height_in);

	// prevent catastrophe
	if( strlen(Align_dir)==0 || strchr(Align_dir,'/')!=NULL || strchr(Align_dir,'.')!=NULL )fatal("align_dir has dangerous value (%s)",Align_dir);
	if( strlen(Correl_dir)==0 || strchr(Correl_dir,'/')!=NULL || strchr(Correl_dir,'.')!=NULL )fatal("correl_dir has dangerous value (%s)",Correl_dir);
}

// remove any previous run and make a new empty directory
static inline void
align_dir_setup()
{
	char	buf[1024];
	int	ret;

	sprintf(buf,"rm -fr %s; mkdir %s",Align_dir,Align_dir);
	ret = system(buf);
	if( ret != 0 )
		fatal("align_dir_setup failed %d",ret);
}

static inline void
correl_dir_setup()
{
	char	buf[1024];
	int	ret;

	sprintf(buf,"rm -fr %s; mkdir %s",Correl_dir,Correl_dir);
	ret = system(buf);
	if( ret != 0 )
		fatal("correl_dir_setup failed %d",ret);
}

static inline void
final_stats()
{
	printf("Aligned %u frames\n",Total_frames);
	printf("Low Trust %u\n",Low_trust);
	printf("Big Dx %u\n",Big_dx);
	printf("Big Dy %u\n",Big_dy);
}

int
main (int argc, char **argv)
{
	setbuf(stdout,NULL);
	configure();
	fft_setup();
	reference_setup();
	check_config();
	align_dir_setup();
	if(Plot_correl)
		correl_dir_setup();
	while(--argc)
		align(*++argv);
	if(Verbose)
		final_stats();
	return 0;
}
