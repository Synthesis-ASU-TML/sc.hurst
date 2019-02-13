//
//  sc.hurst.cpp
//  max-external
//
//  Created by Connor Rawls on 12/17/18.
//

#include "ext.h"                            // standard Max include, always required
#include "ext_obex.h"                       // required for new style Max object
#include "ext_critical.h"                   // for using critical regions
#include <math.h>                           // for log calculations

//#define DEBUG

//=====================OBJECT STRUCT==================
typedef struct _sc_hurst
{
    t_object        ob;
    long series_length;
    long series_max_length;
    long base_division_size;    //the size of the smallest data kernel
    long calc_on_input;
    long show_size_warning;
    //long thread_count;
    double* data_set;
    //t_systhread* threads; //pointer to thread array
#ifdef DEBUG
    long debug;
#endif
    void* out;
    void* out2;
} t_sc_hurst;


//===================HELPER STRUCTS====================

//sent to the helper thread for calculating the mean and standard deviation
typedef struct _sc_hurst_helper_in
{
    double* src_data;
    long idx0;
    long idx1;
    double mean; //filled in by thread (should start as 0)
    double stddev; //filled in by thread (should start as 0)
} t_hurst_helper_ms;

//sent to the helper thread for calculating the range
typedef struct _sc_hurst_helper_range
{
    double* src_data;
    long idx0;
    long idx1;
    double mean;
    double range; //filled in by thread (should start as 0)
} t_hurst_helper_rs;

//sent to linear regression helper function (should be passed as mutable to allow the slope member to be set in the function)
typedef struct _sc_hurst_helper_lin_reg
{
    double* rs; //pointer to array of log2(R(n)/S(n))
    double* size; //pointer to array of log2(block_size)
    long rs_length; //length of rs array
    long size_length; //length of block_size array
    double slope; //slope of the best fit line, filled in by helper function
} t_hurst_helper_lin_reg;

//===================FUNTCTION PROTOTYPES==============

//creation and destruction
void *sc_hurst_new(t_symbol *s, long argc, t_atom *argv);
void sc_hurst_free(t_sc_hurst *x);

//input handling
void sc_hurst_bang(t_sc_hurst *x);
void sc_hurst_int(t_sc_hurst *x, long n);
void sc_hurst_float(t_sc_hurst *x, double f);
void sc_hurst_list(t_sc_hurst *x, t_symbol* a, long argc, t_atom *argv);

//Attribute Mutators
void sc_hurst_set_max_length(t_sc_hurst *x, void *attr, long argc, t_atom *argv);
void sc_hurst_set_length(t_sc_hurst *x, void *attr, long argc, t_atom *argv); //dummy function
void sc_hurst_set_div_size(t_sc_hurst *x, void *attr, long argc, t_atom *argv); //not exposed for now
void sc_hurst_set_coi(t_sc_hurst *x, void *attr, long argc, t_atom *argv);
void sc_hurst_set_size_warning(t_sc_hurst *x, void *attr, long argc, t_atom *argv);
void sc_hurst_set_thread_count(t_sc_hurst *x, void *attr, long argc, t_atom *argv); //not exposed for now

//Attribute Accessors
void sc_hurst_get_max_length(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv);
void sc_hurst_get_length(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv);
void sc_hurst_get_div_size(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv); //not exposed for now
void sc_hurst_get_coi(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv);
void sc_hurst_get_size_warning(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv);
void sc_hurst_get_thread_count(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv); //not exposed for now

//assist function
void sc_hurst_assist(t_sc_hurst *x, void *b, long m, long a, char *s);

//notify for changed attrs from attached objects
t_max_err sc_hurst_notify(t_sc_hurst *x, t_symbol *s, t_symbol *msg, void *sender, void *data);

//general
void sc_hurst_dump(t_sc_hurst *x); //dumps current data set out right outlet
void sc_hurst_get_state(t_sc_hurst *x); //outputs current state of attributes out right outlet
void sc_hurst_clear(t_sc_hurst *x); //clears internal data set

//calculation
void sc_hurst_calculate(t_sc_hurst *x); //calculates hurst exponent if possible

//helper functions for calculations (will only be called by helper threads)
void sc_hurst_stddev_and_mean_helper(t_hurst_helper_ms** x, t_sc_hurst* t); //calculating mean and standard deviation. (requires mutable struct pointer)
void sc_hurst_helper_range(t_hurst_helper_rs** x); //calculating the rescaled range (requires mutable struct pointer)
void sc_hurst_helper_linear_regression(t_hurst_helper_lin_reg** x); //helper function to calculate linear regression of each block size

#ifdef DEBUG
void sc_hurst_set_debug(t_sc_hurst *x, void *attr, long argc, t_atom *argv);
void sc_hurst_get_debug(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv);
#endif

//======================CLASS POINTER VARIABLE============
void *sc_hurst_class;

void ext_main(void *r) {
    t_class *c;
    
    c = class_new("sc.hurst", (method)sc_hurst_new, (method)sc_hurst_free, sizeof(t_sc_hurst), 0L, A_GIMME, 0);
    
    //add functions
    class_addmethod(c, (method)sc_hurst_bang,       "bang",             0);
    class_addmethod(c, (method)sc_hurst_int,        "int",      A_LONG, 0);
    class_addmethod(c, (method)sc_hurst_float,      "float",    A_FLOAT, 0);
    class_addmethod(c, (method)sc_hurst_dump,       "dump",             0);
    class_addmethod(c, (method)sc_hurst_clear,      "clear",            0);
    class_addmethod(c, (method)sc_hurst_get_state,  "getstate",         0);
    class_addmethod(c, (method)sc_hurst_list,       "list",     A_GIMME, 0);
    
    //add attributes
    CLASS_ATTR_LONG(c, "max_length", 0, t_sc_hurst, series_max_length);
    CLASS_ATTR_ACCESSORS(c, "max_length", sc_hurst_get_max_length, sc_hurst_set_max_length);
    
    CLASS_ATTR_LONG(c, "length", ATTR_SET_OPAQUE, t_sc_hurst, series_length);
    CLASS_ATTR_ACCESSORS(c, "length", sc_hurst_get_length, sc_hurst_set_length);
    
    CLASS_ATTR_LONG(c, "calc_on_input", 0, t_sc_hurst, calc_on_input);
    CLASS_ATTR_ACCESSORS(c, "calc_on_input", sc_hurst_get_coi, sc_hurst_set_coi);
    CLASS_ATTR_STYLE(c, "calc_on_input", 0, "onoff");
    
    CLASS_ATTR_LONG(c, "size_warning", 0, t_sc_hurst, show_size_warning);
    CLASS_ATTR_STYLE(c, "size_warning", 0, "onoff");
    CLASS_ATTR_ACCESSORS(c, "size_warning", sc_hurst_get_size_warning, sc_hurst_set_size_warning);
    
#ifdef DEBUG
    CLASS_ATTR_LONG(c, "debug", 0, t_sc_hurst, debug);
    CLASS_ATTR_STYLE(c, "debug", 0, "onoff");
    CLASS_ATTR_ACCESSORS(c, "debug", sc_hurst_get_debug, sc_hurst_set_debug);
#endif
    
    //add assist function
    class_addmethod(c, (method)sc_hurst_assist, "assist", A_CANT, 0);
    
    //add notify function
    class_addmethod(c, (method)sc_hurst_notify, "notify", A_CANT, 0);
    
    //register the module with max
    class_register(CLASS_BOX, c);
    
    sc_hurst_class = c;
}

//Assist function definition
void sc_hurst_assist(t_sc_hurst *x, void *b, long m, long a, char *s) {
    if(m == ASSIST_INLET) {
        sprintf(s, "Floating-point or integer input for hurst exponent calculation");
    } else {
        if(a == 0) {
            sprintf(s, "Outlet %ld : Estimated Hurst Exponent", a);
        } else {
            sprintf(s, "Outlet %ld: dumpout", a);
        }
    }
}

void *sc_hurst_new(t_symbol *s, long argc, t_atom *argv) {
    t_sc_hurst *x = NULL;
    
    if((x = (t_sc_hurst *)object_alloc((t_class*)sc_hurst_class))) {
        //set inital values
        x->series_length = 0;
        x->series_max_length = 256;
        x->base_division_size = 8;
        x->calc_on_input = 1;
        x->show_size_warning = 1;
        x->out = outlet_new(x, 0L);
        x->out2 = outlet_new(x, NULL);
        
        x->data_set = (double*)sysmem_newptr(sizeof(double) * x->series_max_length);
        
        attr_args_process(x, argc, argv);
    } else {
        poststring("Failed to create Hurst");
    }
    
    return (x);
}

//Object destroy function
void sc_hurst_free(t_sc_hurst *x) {
    
    //for now, all we need to do is get rid of the data set pointers.
    if(x->data_set != NULL) { //don't try to free non-existent data
        double* temp = x->data_set;
        for(int i = 0; (i < x->series_length || i < x->series_max_length); i++) {
            double* t2 = temp;
            temp++; //advance before we free the pointer we're referencing
            sysmem_freeptr(t2);
        }
    }
}

//notify for changed attrs from attached objects
t_max_err sc_hurst_notify(t_sc_hurst *x, t_symbol *s, t_symbol *msg, void *sender, void *data) {
    //this is just defined for safety at the moment
    return 0;
}

/*=============================================================
 ==================INPUT HANDLERS==============================
 ==============================================================*/

void sc_hurst_bang(t_sc_hurst *x){ //just attempt to calculate on bang
    sc_hurst_calculate(x);
}

void sc_hurst_int(t_sc_hurst *x, long n){ //add data to the array

    critical_enter(0); //block until done editing array
    if(x->series_length == x->series_max_length) {
        double* temp = x->data_set; //temporary pointer to data set
        double* temp2 = temp; //pointer to index 1 of set
        temp2++;
        
        sysmem_copyptr(temp2, temp, sizeof(double) * (x->series_max_length - 1)); //shift data 1 index to the left
        //for(int i = 0; i < (x->series_max_length - 1); i++, temp++){} //move to the end of the array
        temp += x->series_max_length - 1;
        *temp = (double)n; //fill last index with new data
        
    } else {
        double* temp = x->data_set + x->series_length; //get a temporary pointer to the data set
        //for(int i = 0; i < x->series_length; i++) {temp++;} //increment to next available index
        *temp = (double)n; //add data
        x->series_length++;
    }
    critical_exit(0);
    
    if(x->calc_on_input == 1) {
        sc_hurst_calculate(x);
    }
}

void sc_hurst_float(t_sc_hurst *x, double f) {
    critical_enter(0); //block until done editing array
    if(x->series_length == x->series_max_length) {
        double* temp = x->data_set; //temporary pointer to data set
        double* temp2 = temp; //pointer to index 1 of set
        temp2++;
        
        sysmem_copyptr(temp2, temp, sizeof(double) * (x->series_max_length - 1)); //shift data 1 index to the left
        //for(int i = 0; i < (x->series_max_length - 1); i++, temp++){} //move to the end of the array
        temp += x->series_max_length - 1;
        *temp = f; //fill last index with new data
        
    } else {
        double* temp = x->data_set + x->series_length; //get a temporary pointer to the data set
        //for(int i = 0; i < x->series_length; i++) {temp++;} //increment to next available index
        *temp = f; //add data
        x->series_length++;
    }
    critical_exit(0);
    
    if(x->calc_on_input == 1) {
        sc_hurst_calculate(x);
    }
}

void sc_hurst_list(t_sc_hurst *x, t_symbol* a, long argc, t_atom *argv) {
    t_atom* arg_temp = argv;
    for(int i = 0; i < argc; i++, arg_temp++) {
        switch(atom_gettype(arg_temp)) {
            case A_LONG:
                break;
            case A_FLOAT:
                break;
            default:
                object_warn((t_object*)x, "Received non-numeric input");
                return;
        }
    }
    
    arg_temp = argv;
    int idx = 0;
    double* data_list;
    long data_size = argc;
    if(argc > x->series_max_length) {
        arg_temp += argc - x->series_max_length;
        idx = argc - x->series_max_length;
        data_size = x->series_max_length;
    }
    
    data_list = (double*)sysmem_newptr(sizeof(double) * data_size);
    
    double* data_temp = data_list;
    
    for(; idx < argc; idx++, arg_temp++, data_temp++) {
        switch(atom_gettype(arg_temp)) {
            case A_LONG:
                *data_temp = (double)atom_getlong(arg_temp);
                break;
            case A_FLOAT:
                *data_temp = atom_getfloat(arg_temp);
                break;
        }
    }
    
    long tot_size = x->series_length + data_size;
    
    long del_idx = 0;
    
    if(tot_size > x->series_max_length) {
        del_idx = tot_size - x->series_max_length;
        double* temp0 = x->data_set + del_idx;
        double* temp1 = x->data_set;
        sysmem_copyptr(temp0, temp1, sizeof(double) * (x->series_max_length - del_idx));
        x->series_length = (x->series_length - del_idx > 0) ? (x->series_length - del_idx) : 0;
    }
    
    double* temp = x->data_set + x->series_length;
    
    sysmem_copyptr(data_list, temp, sizeof(double) * data_size);
    
    //free data
    data_temp = data_list;
    for(int i = 0; i < data_size; i++) {
        double* d2 = data_temp;
        data_temp++;
        sysmem_freeptr(d2);
    }
    
    x->series_length += data_size;
}

/*==================================================================================
 ========================ATTRIBUTE MUTATORS=========================================
 ====================================================================================*/

void sc_hurst_set_max_length(t_sc_hurst *x, void *attr, long argc, t_atom *argv){
    if(argc && argv) {
        long temp_sl = 0;
        
        switch(atom_gettype(argv)){
            case A_LONG:
                temp_sl = atom_getlong(argv);
                break;
            case A_FLOAT:
                temp_sl = (long)atom_getfloat(argv);
                break;
            default:
                object_error((t_object *)x, "Bad value for max_length. Expected a positive integer");
                return;
                break;
        }
        
        if(temp_sl > 16) {
            if(temp_sl < x->series_max_length) {
                double* temp = (double*)sysmem_newptr(sizeof(double) * temp_sl);
                double* t2 = x->data_set; //make sure to only copy the newest data
                if(temp_sl < x->series_length) {
                    t2 += x->series_length - temp_sl;
                }
                sysmem_copyptr(t2, temp, sizeof(double) * temp_sl);
                
                double* d1 = x->data_set; //create pointer to data we will free
                x->data_set = temp; //point struct pointer to newly copied data
                
                //free pointers
                for(int i = 0; i < x->series_max_length; i++) {
                    double* d2 = d1;
                    d1++;
                    sysmem_freeptr(d2);
                }
                
                x->series_max_length = temp_sl;
                
            } else if(temp_sl > x->series_max_length) {
                double* temp = (double*)sysmem_newptr(sizeof(double) * temp_sl); //allocate memory of the new size
                sysmem_copyptr(x->data_set, temp, sizeof(double) * x->series_max_length); //copy existing data
                double* d1 = x->data_set; //make pointer for freeing old data
                x->data_set = temp; //struct pointer to new data
                
                //free pointers
                for(int i = 0; i < x->series_max_length; i++) {
                    double* d2 = d1;
                    d1++;
                    sysmem_freeptr(d2);
                }
                x->series_max_length = temp_sl;
                
            } else {
                //fail silently and do nothing
            }
        } else {
            object_error((t_object *)x, "Bad value for max_length. Expected a poisitve integer >= 16");
        }
    }
}

void sc_hurst_set_length(t_sc_hurst *x, void *attr, long argc, t_atom *argv) {
    /*
     I'm here, and you can call me, but I don't do anything.
     */
}
void sc_hurst_set_div_size(t_sc_hurst *x, void *attr, long argc, t_atom *argv) {
    /*
     I'll eventually do something, but for now I do nothing
     */
}
void sc_hurst_set_thread_count(t_sc_hurst *x, void *attr, long argc, t_atom *argv){
    /*
     One day I'll become a real boy
     */
}

void sc_hurst_set_coi(t_sc_hurst *x, void *attr, long argc, t_atom *argv){
    if(argc && argv) {
        long temp_coi = atom_getlong(argv);
        
        switch(atom_gettype(argv)) {
            case A_LONG:
                temp_coi = atom_getlong(argv);
                break;
            case A_FLOAT:
                temp_coi = (long)atom_getfloat(argv);
                break;
            default:
                object_error((t_object *)x, "bad value received for calc_on_input");
                return;
                break;
        }
        
        if(temp_coi > 1){temp_coi = 1;}
        if(temp_coi < 0){temp_coi = 0;}
        x->calc_on_input = temp_coi;
    }
}

void sc_hurst_set_size_warning(t_sc_hurst *x, void *attr, long argc, t_atom *argv) {
    if(argc && argv) {
        long temp_sw = 0;
        
        switch(atom_gettype(argv)) {
            case A_LONG:
                temp_sw = atom_getlong(argv);
                break;
            case A_FLOAT:
                temp_sw = (long)atom_getfloat(argv);
                break;
            default:
                object_error((t_object *)x, "bad value received for size_warning");
                return;
                break;
        }
        if(temp_sw > 1) {temp_sw = 1;}
        if(temp_sw < 0) {temp_sw = 0;}
        
        x->show_size_warning = temp_sw;
    }
}

/*==================================================================================
 ========================ATTRIBUTE ACCESSORS=========================================
 ====================================================================================*/

void sc_hurst_get_max_length(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv){
    char alloc;
    long sl = 0;
    
    atom_alloc(argc, argv, &alloc);
    sl = x->series_max_length;
    atom_setlong(*argv, sl);
}

void sc_hurst_get_length(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv){
    char alloc;
    long csize = 0;
    
    atom_alloc(argc, argv, &alloc);
    csize = x->series_length;
    atom_setlong(*argv, csize);
}

void sc_hurst_get_div_size(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv){
    /*
     One day I'll be a realy boy, for now I'm just a copy of the function above me.
     */
    char alloc;
    long csize = 8;
    
    atom_alloc(argc, argv, &alloc);
    //csize = x->series_length;
    atom_setlong(*argv, csize);
}

void sc_hurst_get_coi(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv){
    char alloc;
    long coi = 0;
    
    atom_alloc(argc, argv, &alloc);
    coi = x->calc_on_input;
    atom_setlong(*argv, coi);
}

void sc_hurst_get_size_warning(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv) {
    char alloc;
    long sw = 0;
    
    atom_alloc(argc, argv, &alloc);
    sw = x->show_size_warning;
    atom_setlong(*argv, sw);
}

void sc_hurst_get_thread_count(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv) {
    /*
     One day I'll be a real boy, for now I'm just a copy of the function above me.
     */
    char alloc;
    long sw = 1;
    
    atom_alloc(argc, argv, &alloc);
    //sw = x->show_size_warning;
    atom_setlong(*argv, sw);
}

/*==================================================================================
 ========================GENERAL FUNCTIONS=========================================
 ====================================================================================*/

void sc_hurst_dump(t_sc_hurst *x) {
    critical_tryenter(0);
    
    if(x->series_length > 0){
        double* d = x->data_set;
        
        void* mem = sysmem_newptr(sizeof(t_atom) * (x->series_length + 1));
        t_atom* list = (t_atom*)mem;
        t_atom* temp_list = list;
        atom_setsym(temp_list, gensym("values"));
        temp_list++;
        for(int i = 0; i < x->series_length; i++, d++, temp_list++) {
            atom_setfloat(temp_list, *d);
        }
        outlet_list((void*)x->out, gensym("values"), x->series_length, list);
        
        sysmem_freeptr(mem);
        
    }
    critical_exit(0);
}

void sc_hurst_get_state(t_sc_hurst *x){
    
    void* state = sysmem_newptr(sizeof(t_atom) * 2);
    
    //calc on input
    t_atom* coi_list = (t_atom*)state;
    atom_setsym(coi_list, gensym("calc_on_input"));
    coi_list++;
    atom_setlong(coi_list, x->calc_on_input);
    outlet_list(x->out, gensym("calc_on_input"), 2, (t_atom*)state);
    coi_list = NULL;
    
    //max series length
    t_atom* temp_list = (t_atom*)state;
    atom_setsym(temp_list, gensym("max_length"));
    temp_list++;
    atom_setlong(temp_list, x->series_max_length);
    outlet_list(x->out, gensym("max_length"), 2, (t_atom*)state);
    
    //current series length
    temp_list = (t_atom*)state;
    atom_setsym(temp_list, gensym("length"));
    temp_list++;
    atom_setlong(temp_list, x->series_length);
    outlet_list(x->out, gensym("length"), 2, (t_atom*)state);
    
    //size warning
    temp_list = (t_atom*)state;
    atom_setsym(temp_list, gensym("size_warning"));
    temp_list++;
    atom_setlong(temp_list, x->show_size_warning);
    outlet_list(x->out, gensym("size_warning"), 2, (t_atom*)state);
    
    temp_list = NULL;
    sysmem_freeptr(state);
    
    sc_hurst_dump(x);
}

void sc_hurst_clear(t_sc_hurst *x){
    critical_tryenter(0);
    double* emptyd = (double*)sysmem_newptr(sizeof(double) * x->series_max_length);
    sysmem_copyptr(emptyd, x->data_set, sizeof(double) * x->series_max_length);
    for(int i = 0; i < x->series_max_length; i++) {
        double* temp = emptyd;
        emptyd++;
        sysmem_freeptr(temp);
    }
    x->series_length = 0;
    critical_exit(0);
}


/*==================================================================================
 ========================CALCULATION FUNCTIONS=========================================
 ====================================================================================*/

void sc_hurst_calculate(t_sc_hurst *x) {
    //escape if there are not enough values to calculate
    if(x->series_length < 16) {
        if(x->show_size_warning == 1) {
            object_warn((t_object*)x, "Too few values to calculate Hurst Exponent.");
            object_warn((t_object*)x, "Requires 16 values, currently have %ld.", x->series_length);
        }
        return;
    }
    
    //object_post((t_object*)x, "Beginning Calculations");
    
    long div_size = 2;
    
    if(x->series_length < 64) {
        div_size = 2;
    } else if(x->series_length < 128) {
        div_size = 4;
    } else if(x->series_length < 256) {
        div_size = 6;
    } else {
        div_size = 8;
    }
    
    //object_post((t_object*)x, "Base division size: %ld", div_size);
    
    long layer_count = 0;
    for(int i = 0; pow(2, i) * div_size <= x->series_length ; i++) {
        //long test_val = pow(2, i) * div_size;
        //object_post((t_object *)x, "Current step: %ld", test_val);
        layer_count++;
    }
    
   //object_post((t_object*)x, "Layer count: %ld", layer_count);
    
    double* rs_avg = (double*)sysmem_newptr(sizeof(double) * layer_count);
    double* size = (double*)sysmem_newptr(sizeof(double) * layer_count);
    
    double* rs_temp = rs_avg;
    double* size_temp = size;
    
    for(int i = 0; i < layer_count; i++) {
        long cur_div_size = pow(2, i) * div_size;
        
        long layer_size = x->series_length / cur_div_size;
        
        //object_post((t_object*)x, "Layer Block Count: %ld", layer_size);
        //object_post((t_object*)x, "Layer Div Size: %ld", cur_div_size);
        double rs_sum = 0;
        for(int j = 0; j < layer_size; j++) {
            //create and fill helper struct pointer
            t_hurst_helper_ms* ms_temp = (t_hurst_helper_ms*)sysmem_newptr(sizeof(t_hurst_helper_ms));
            ms_temp->src_data = x->data_set;
            ms_temp->idx0 = (j * (pow(2, i) * div_size));
            long end_idx = ((j + 1) * (pow(2, i) * div_size));
            ms_temp->idx1 = (end_idx < x->series_length) ? end_idx : (x->series_length - 1);
            ms_temp->mean = 0; //filled in function
            ms_temp->stddev = 0; //filled in function
            
            sc_hurst_stddev_and_mean_helper(&ms_temp, x); //compute mean and standard deviation for the block.
            
            //copy out standard deviation and mean
            double stddev = (ms_temp->stddev > 0) ? ms_temp->stddev : 0.0001;
            double mean = ms_temp->mean;
            //free the helper pointer after removing the reference to the struct data
            ms_temp->src_data = NULL;
            sysmem_freeptr(ms_temp);

#ifdef DEBUG
            if(x->debug > 0) {
                object_post((t_object*)x, "div size: %ld, block: %ld, mean: %f, stddev: %f", cur_div_size, j, mean, stddev);
            }
#endif
            
            t_hurst_helper_rs* rsa_temp = (t_hurst_helper_rs*)sysmem_newptr(sizeof(t_hurst_helper_rs));
            rsa_temp->src_data = x->data_set;
            rsa_temp->mean = mean;
            rsa_temp->idx0 = (j * (pow(2, i) * div_size));
            rsa_temp->idx1 = (end_idx < x->series_length) ? end_idx : (x->series_length - 1);
            rsa_temp->range = 0;
            
            sc_hurst_helper_range(&rsa_temp);
            double range = rsa_temp->range;
            double rs = range / stddev;
            rs_sum += rs;
            
#ifdef DEBUG
            if(x->debug > 0) {
                object_post((t_object*)x, "range: %f, rs: %f", range, rs);
            }
#endif
            
            //free the struct pointer
            rsa_temp->src_data = NULL;
            sysmem_freeptr(rsa_temp);
        }
        
        *rs_temp = log2(rs_sum / ((layer_size > 0) ? layer_size : 1));

        *size_temp = log2((pow(2, i) * div_size));

        
#ifdef DEBUG
        if(x->debug > 0) {
            object_post((t_object*)x, "log2rs: %f, log2size: %f", *rs_temp, *size_temp);
        }
#endif
        rs_temp++;
        size_temp++;
    }
    
    //create struct for linear regression helper function
    t_hurst_helper_lin_reg* lr_temp = (t_hurst_helper_lin_reg*)sysmem_newptr(sizeof(t_hurst_helper_lin_reg));
    lr_temp->rs = rs_avg;
    lr_temp->size = size;
    lr_temp->rs_length = layer_count;
    lr_temp->size_length = layer_count;
    lr_temp->slope = 0;
    
    sc_hurst_helper_linear_regression(&lr_temp);
    
    double hurst_exp = lr_temp->slope;
    
    //cleanup allocated memory
    lr_temp->rs = NULL;
    lr_temp->size = NULL;
    sysmem_freeptr(lr_temp);
    
    for(int i = 0; i < layer_count; i++){
        double* r2 = rs_avg;
        double* s2 = size;
        rs_avg++;
        size++;
        sysmem_freeptr(r2);
        sysmem_freeptr(s2);
    }
    
    //output the computed data
    outlet_float(x->out2, hurst_exp);
}


void sc_hurst_stddev_and_mean_helper(t_hurst_helper_ms** x, t_sc_hurst* t){
    double* temp = (*x)->src_data + (*x)->idx0;
    //for(int i = 0; i < (*x)->idx0; i++, temp++) {} //move to idx0
    double* t2 = temp;
    
    long length = (*x)->idx1 - (*x)->idx0;
    
    //calculate mean
    double sum = 0;
    for(int i = 0; i < length; i++, t2++) {
        sum += *t2;
    }
    
#ifdef DEBUG
    if(t->debug > 0) {
        object_post((t_object*)t, "length: %ld, sum: %f", length, sum);
    }
#endif
    
    double mean = sum / ((double)length);
    
    //calculate standard deviation
    double stddev_sum = 0;
    for(int i = 0; i < length; i++, temp++) {
        stddev_sum += pow((*temp - mean), 2);
    }
    
    double stddev = pow((stddev_sum * (1 / ((double)length))), 0.5);
    
    (*x)->mean = mean;
    (*x)->stddev = stddev;
    
#ifdef DEBUG
    if(t->debug > 0) {
        object_post((t_object*)t, "idx0: %ld, idx1: %ld, mean: %f, stddev: %f", (*x)->idx0, (*x)->idx1, mean, stddev);
    }
#endif
}

void sc_hurst_helper_range(t_hurst_helper_rs** x) {
    double* temp = (*x)->src_data + (*x)->idx0;
    //for(int i = 0; i < (*x)->idx0; i++, temp++){} //move to the starting index
    long length = (*x)->idx1 - (*x)->idx0;
    
    double min = *temp - (*x)->mean;
    double max = *temp - (*x)->mean;
    
    double t = 0;
    for(int i = 0; i < length; i++, temp++) {
        t += (*temp - (*x)->mean); //t'(n) = (t(n) - mean) + t'(n-1)
        if(t > max) {
            max = t;
        } else if(t < min) {
            min = t;
        }
    }
    (*x)->range = max - min;
}

void sc_hurst_helper_linear_regression(t_hurst_helper_lin_reg** x){
    double* s = (*x)->size; //x-axis values
    double* r = (*x)->rs; //y-axis values
    double length = (*x)->rs_length;
    double sumx = 0;
    double sumy = 0;
    double sumxsq = 0;
    double sumxy = 0;
    
    for(int i = 0; i < length; i++, s++, r++){
        sumx += *s;
        sumy += *r;
        sumxsq += pow(*s, 2);
        sumxy += *s * *r;
    }
    
    double denom = (length * sumxsq) - pow(sumx, 2);
    double slope = ((length * sumxy) - (sumx * sumy)) / denom;
    
    (*x)->slope = slope;
}

/*================================================================
 =========================DEBUGGING ATTRIBUTE=====================
 =================================================================*/

#ifdef DEBUG
void sc_hurst_set_debug(t_sc_hurst *x, void *attr, long argc, t_atom *argv){
    if(argc && argv) {
        long temp_d = atom_getlong(argv);
        
        switch(atom_gettype(argv)) {
            case A_LONG:
                temp_d = atom_getlong(argv);
                break;
            case A_FLOAT:
                temp_d = (long)atom_getfloat(argv);
                break;
            default:
                object_error((t_object *)x, "bad value received for debug");
                return;
                break;
        }
        
        if(temp_d > 1){temp_d = 1;}
        if(temp_d < 0){temp_d = 0;}
        x->debug = temp_d;
    }
}
void sc_hurst_get_debug(t_sc_hurst *x, t_object *attr, long *argc, t_atom **argv){
    char alloc;
    long d = 0;
    
    atom_alloc(argc, argv, &alloc);
    d = x->debug;
    atom_setlong(*argv, d);
}
#endif
