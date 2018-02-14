/*
	@file
	timestruct - extracts the amplitude timing structure from an audio buffer

    todo:
    *create a copy of buffer to analyze and calculate on it, so i dont have to modify the original
    *handle stereo buffers
    *smooth out jitter in peak/dip detection
*/

#include "ext.h"
#include "z_dsp.h"
#include "math.h"
#include "ext_buffer.h"
#include "ext_atomic.h"
#include "ext_obex.h"


//*************************************************************************************************//
//the objects datastructure
typedef struct _timestruct
{
    t_object w_obj; //the structure for the head of any object which wants to have inlets or outlets
	t_buffer_ref *w_buf; //reference to a buffer~ object in Max
	t_symbol *w_name; //the name of the buffer~ object
	long w_len; //number of frames
    double w_msr; //global sample rate
    t_atom w_collList[7];
    t_atom w_scoreList[19];
	t_bool w_buffer_modified;
    void *w_bangOutlet;
    void *w_collOutlet;
    void *w_scoreOutlet;
} t_timestruct;


//*************************************************************************************************//
//function prototypes
void *timestruct_new(t_symbol *s,  long argc, t_atom *argv);
void timestruct_free(t_timestruct *x);
t_max_err timestruct_notify(t_timestruct *x, t_symbol *s, t_symbol *msg, void *sender, void *data);
void timestruct_assist(t_timestruct *x, void *b, long m, long a, char *s);
void timestruct_set(t_timestruct *x, t_symbol *s, long ac, t_atom *av);
void timestruct_bang(t_timestruct *x);
void timestruct_float(t_timestruct *x, double f);
void timestruct_int(t_timestruct *x, long n);
void timestruct_dblclick(t_timestruct *x);

void timestruct_getTime(t_timestruct *x, t_symbol *s, long ac, t_atom *av);
void timestruct_doGetTime(t_timestruct *x, double factor);
void timestruct_sendOutlets(t_timestruct *x, double *tripletList, long collIndex);


static t_symbol *ps_buffer_modified;
static t_class *s_timestruct_class;


//*************************************************************************************************//
void ext_main(void *r)
{
	t_class *c = class_new("timestruct", (method)timestruct_new, (method)timestruct_free, sizeof(t_timestruct), NULL, A_GIMME, 0);

    class_addmethod(c, (method)timestruct_bang, "bang", 0);
	class_addmethod(c, (method)timestruct_float, "float", A_FLOAT, 0);
	class_addmethod(c, (method)timestruct_int, "int", A_LONG, 0);
	class_addmethod(c, (method)timestruct_set, "set", A_GIMME, 0);
	class_addmethod(c, (method)timestruct_getTime, "getTime", A_GIMME, 0);
	class_addmethod(c, (method)timestruct_assist, "assist", A_CANT, 0);
	class_addmethod(c, (method)timestruct_dblclick, "dblclick", A_CANT, 0);
	class_addmethod(c, (method)timestruct_notify, "notify", A_CANT, 0);
	class_addmethod(c, (method)timestruct_doGetTime, "doGet", A_CANT, 0);
    class_addmethod(c, (method)timestruct_sendOutlets, "sendOutlets", A_CANT, 0);
    
	class_register(CLASS_BOX, c);
	s_timestruct_class = c;

	ps_buffer_modified = gensym("buffer_modified");
}


//*************************************************************************************************//
void *timestruct_new(t_symbol *s,  long argc, t_atom *argv)
{
	t_timestruct *x = (t_timestruct *)object_alloc(s_timestruct_class);
	t_symbol *buf=0;
	x->w_msr = sys_getsr(); //Query MSP for the global sample rate.
    
	buf = atom_getsymarg(0,argc,argv);

	x->w_name = buf;
    x->w_bangOutlet = bangout((t_object *)x);
    x->w_collOutlet = listout((t_object *)x);
    x->w_scoreOutlet = outlet_new((t_object *)x, NULL);
    
	//create a new buffer reference, initially referencing a buffer with the provided name
	x->w_buf = buffer_ref_new((t_object *)x, x->w_name);
    //store number of frames of buffer
    x->w_len = buffer_getframecount(buffer_ref_getobject(x->w_buf));

    object_post((t_object *)x, "timestruct.mxo -- lukas nowok -- v0.1 - 26.Jan.2017");
	return (x);
}


//*************************************************************************************************//
void timestruct_free(t_timestruct *x)
{
	//must free buffer reference when no longer used
	object_free(x->w_buf);
}


//*************************************************************************************************//
void timestruct_getTime(t_timestruct *x, t_symbol *s, long ac, t_atom *av)
{
    timestruct_doGetTime(x, atom_getfloat(av));
}


//*************************************************************************************************//
void timestruct_doGetTime(t_timestruct *x, double factor)
{
    t_float *tab; //pointer to the first sample in memory
    t_buffer_obj *buffer = buffer_ref_getobject(x->w_buf);
    
    //claim the buffer∼ and get a pointer to the first sample in memory
    tab = buffer_locksamples(buffer);
    
    t_float maxSample = 0.;
    t_float normalizeFactor = 1.;

    //2 x expontial smoothing for a smoooooth amplitude contour, yum
    //t_float factor = 0.001; //smoothingfactor
    t_float average; //working average
    for(int z = 0; z < 2; z++)
    {
        average = fabsf(tab[0]);
        for(int i = 1; i < x->w_len; i++)
        {
            //workingAverage = (newValue*smoothingFactor) + ( ( 1.0 - smoothingFactor) * workingAverage)
            average = fabsf(tab[i])*factor + (1-factor)*average;
            //write the average (smoothed) value back into the buffer
            tab[i] = average;
            
            //find maximum sample peak
            if(maxSample < tab[i])
            {
                maxSample = tab[i];
                normalizeFactor = 1/maxSample;
            };
        };
    };
    
    //normalizing
    for(int i = 0; i < x->w_len; i++)
    {
        tab[i] = tab[i]*normalizeFactor;
    };

    //
    //detect ramp up/ramp down phases (peaks/dips)
    //
    t_int distance = 1000;
    t_float first = tab[0];
    t_float last = tab[distance];
    t_int collIndex = 0;
    double tripletList[6]; //[posDip, ampDip, posPeak, ampPeak, posDip, ampDip]
    t_int tripletIndex = 0;
    bool up = (first < last);
    
    //clear the bach.roll before writing into it
    outlet_anything(x->w_scoreOutlet, gensym("clear"),0,NIL);
    
    for(int i = 0; i < x->w_len-distance; i++)
    {
        first = tab[i];
        last = tab[i+distance];
        if(first < last && up)
        {//dip
            if(tripletIndex==0)
            {
                tripletList[0] = i;
                tripletList[1] = first;
                tripletIndex = (tripletIndex+1)%3;
            }
            else if(tripletIndex==2)
            {
                tripletList[4] = i;
                tripletList[5] = first;
                tripletIndex = (tripletIndex+1)%3;
                
                //call timestruct_sendOutlets
                timestruct_sendOutlets(x, tripletList, collIndex);
                collIndex += 1; //count up
                
                //last dip of previous is fist of... next, burp!
                tripletList[0] = i;
                tripletList[1] = first;
                tripletIndex = (tripletIndex+1)%3;
            };
            //timestruct_sendOutlets(x, 0, i, collIndex, first);
            up = !up;
        }
        else if(first > last && !up)
        {//peak
            //if there is a dip before this peak, store its data (pos and amp)
            if(tripletIndex==1)
            {
                tripletList[2] = i-distance;
                tripletList[3] = first;
                tripletIndex = (tripletIndex+1)%3;
            }
            //timestruct_sendOutlets(x, 1, i, collIndex, first);
            up = !up;
        };
    };
    
    //release claim on the buffer∼ contents so that other objects may read/write to the buffer∼
    buffer_unlocksamples(buffer);
    //mark buffer as dirty so that objects such as waveform∼ update their rendering of the contents
    object_method(buffer, gensym("dirty"));
    
    //bang rightmost outlet when done
    outlet_bang(x->w_bangOutlet);
}


//*************************************************************************************************//
void timestruct_sendOutlets(t_timestruct *x, double *tripletList, long collIndex)
{
    //tripletList = [posDip, ampDip, posPeak, ampPeak, posDip, ampDip]
    
    //send output to collOutlet
    atom_setlong(x->w_collList, collIndex); //index into the coll list
    for(int i=0; i<3; i++)
    {
        atom_setlong(x->w_collList+((i*2)+1), (tripletList[i*2]/x->w_msr)*1000); //position as offset in ms
        atom_setfloat(x->w_collList+((i*2)+2), tripletList[(i*2)+1]); //amplitude of the peak/dip
    };
    outlet_list(x->w_collOutlet,0L,7,x->w_collList);

    
    //send output to scoreOutlet
    //format: addchord(offset(pitch duration velocity (slots (1 (x y slope)(x y slope)(x y slope)))))
    atom_setsym(x->w_scoreList, gensym("("));
    atom_setlong(x->w_scoreList+1, (tripletList[0]/x->w_msr)*1000); //offset in ms of the dip ('start of the event...')
    atom_setsym(x->w_scoreList+2, gensym("("));
    atom_setlong(x->w_scoreList+3, 6000); //pitch
    atom_setlong(x->w_scoreList+4, ((tripletList[4]-tripletList[0])/x->w_msr)*1000); //duration (might be problematic for quantisation)
    //atom_setlong(x->w_scoreList+4, 100); //duration
    atom_setlong(x->w_scoreList+5, tripletList[3]*128); //velocity
    
    //////////////////
    //slot information
    atom_setsym(x->w_scoreList+6, gensym("(slots(1("));
    //slot1, point 1
    atom_setlong(x->w_scoreList+7, 0); //position
    atom_setfloat(x->w_scoreList+8, tripletList[1]); //amplitude
    atom_setlong(x->w_scoreList+9, 0); //slope
    atom_setsym(x->w_scoreList+10, gensym(")("));
    
    //slot1, point 2
    atom_setfloat(x->w_scoreList+11, (tripletList[2]-tripletList[0])/(tripletList[4]-tripletList[0])); //peak position
    atom_setfloat(x->w_scoreList+12, tripletList[3]); //amplitude
    atom_setlong(x->w_scoreList+13, 0); //slope
    atom_setsym(x->w_scoreList+14, gensym(")("));
    
    //slot1, point 3
    atom_setlong(x->w_scoreList+15, 1); //position
    atom_setfloat(x->w_scoreList+16, tripletList[5]); //amplitude
    atom_setlong(x->w_scoreList+17, 0); //slope
    atom_setsym(x->w_scoreList+18, gensym(")))))"));
    
    //send it all out
    outlet_anything(x->w_scoreOutlet, gensym("addchord"), 19, x->w_scoreList);
}


//*************************************************************************************************//
//a notify method is required for the buffer reference
//this handles notifications when the buffer appears, disappears, or is modified.
t_max_err timestruct_notify(t_timestruct *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
	if (msg == ps_buffer_modified)
		x->w_buffer_modified = true;
    
    //update numFrames
    x->w_len = buffer_getframecount(buffer_ref_getobject(x->w_buf));
    //notify method for the #t_buffer_ref
	return buffer_ref_notify(x->w_buf, s, msg, sender, data);
}


//*************************************************************************************************//
void timestruct_doset(t_timestruct *x, t_symbol *s, long ac, t_atom *av)
{
	t_symbol *name;

	name = (ac) ? atom_getsym(av) : gensym("");

    //change the buffer used by buffer reference
	buffer_ref_set(x->w_buf, name);
    //update numFrames
    x->w_len = buffer_getframecount(buffer_ref_getobject(x->w_buf));
}


//*************************************************************************************************//
//calls to set the buffer ref should happen on the main thread only
void timestruct_set(t_timestruct *x, t_symbol *s, long ac, t_atom *av)
{
	defer(x, (method)timestruct_doset, s, ac, av);
}


//*************************************************************************************************//
void timestruct_float(t_timestruct *x, double f)
{
    //timestruct_getTiming(x,f);
}


//*************************************************************************************************//
void timestruct_bang(t_timestruct *x)
{
    //timestruct_getTiming(x);
}


//*************************************************************************************************//
void timestruct_int(t_timestruct *x, long n)
{
}


//*************************************************************************************************//
void timestruct_dblclick(t_timestruct *x)
{
	buffer_view(buffer_ref_getobject(x->w_buf));
}


//*************************************************************************************************//
//messages when mousing over in- and outlets
void timestruct_assist(t_timestruct *x, void *b, long m, long a, char *s)
{
    if (m == ASSIST_INLET)
    {	//inlets
        switch (a)
        {
            case 0:	snprintf_zero(s, 256, "bang to output time structure"); break;
        }
    };
    
    if (m == ASSIST_OUTLET)
    {	//outlet
        switch (a)
        {
            case 0: snprintf_zero(s, 256, "(score data) bach.roll format output"); break;
            case 1: snprintf_zero(s, 256, "(raw time data) coll format output"); break;
            case 2: snprintf_zero(s, 256, "bang when done reading"); break;
        }
    }
}








