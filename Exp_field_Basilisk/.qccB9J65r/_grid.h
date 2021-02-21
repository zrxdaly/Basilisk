size_t datasize = 16*sizeof (double);
static int metric (const int i, const double t, Event * _ev);
static int metric_expr0 (int * ip, double * tp, Event * _ev);
static int defaults (const int i, const double t, Event * _ev);
static int defaults_expr0 (int * ip, double * tp, Event * _ev);
static int init (const int i, const double t, Event * _ev);
static int init_expr0 (int * ip, double * tp, Event * _ev);
static int set_dtmax (const int i, const double t, Event * _ev);
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev);
static int stability (const int i, const double t, Event * _ev);
static int stability_expr0 (int * ip, double * tp, Event * _ev);
static int vof (const int i, const double t, Event * _ev);
static int vof_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection (const int i, const double t, Event * _ev);
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_diffusion (const int i, const double t, Event * _ev);
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev);
static int properties (const int i, const double t, Event * _ev);
static int properties_expr0 (int * ip, double * tp, Event * _ev);
static int advection_term (const int i, const double t, Event * _ev);
static int advection_term_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term (const int i, const double t, Event * _ev);
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration (const int i, const double t, Event * _ev);
static int acceleration_expr0 (int * ip, double * tp, Event * _ev);
static int projection (const int i, const double t, Event * _ev);
static int projection_expr0 (int * ip, double * tp, Event * _ev);
static int end_timestep (const int i, const double t, Event * _ev);
static int end_timestep_expr0 (int * ip, double * tp, Event * _ev);
static int adapt (const int i, const double t, Event * _ev);
static int adapt_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_0 (const int i, const double t, Event * _ev);
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection_0 (const int i, const double t, Event * _ev);
static int tracer_advection_0_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_diffusion_0 (const int i, const double t, Event * _ev);
static int tracer_diffusion_0_expr0 (int * ip, double * tp, Event * _ev);
static int properties_0 (const int i, const double t, Event * _ev);
static int properties_0_expr0 (int * ip, double * tp, Event * _ev);
static int init_0 (const int i, const double t, Event * _ev);
static int init_0_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_0 (const int i, const double t, Event * _ev);
static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_diffusion_1 (const int i, const double t, Event * _ev);
static int tracer_diffusion_1_expr0 (int * ip, double * tp, Event * _ev);
static int mov (const int i, const double t, Event * _ev);
static int mov_expr0 (int * ip, double * tp, Event * _ev);
static int adapt_0 (const int i, const double t, Event * _ev);
static int adapt_0_expr0 (int * ip, double * tp, Event * _ev);
static int end (const int i, const double t, Event * _ev);
static int end_expr0 (int * ip, double * tp, Event * _ev);
static void _set_boundary0 (void);
static void _set_boundary1 (void);
static void _set_boundary2 (void);
static void _set_boundary3 (void);
static void _set_boundary4 (void);
static void _set_boundary5 (void);
static void _set_boundary6 (void);
static void _set_boundary7 (void);
static void _set_boundary8 (void);
static void _set_boundary9 (void);
static void _set_boundary10 (void);
static void _set_boundary11 (void);
void _init_solver (void) {
  void init_solver();
  init_solver();
  Events = (Event *) pmalloc (sizeof (Event), __func__, __FILE__, __LINE__);
  Events[0].last = 1;
  event_register ((Event){ 0, 1, metric, {metric_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/embed.h", 891, "metric"});
  event_register ((Event){ 0, 1, defaults, {defaults_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 126, "defaults"});
  event_register ((Event){ 0, 1, defaults_0, {defaults_0_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/tracer.h", 25, "defaults"});
  event_register ((Event){ 0, 1, init, {init_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 181, "init"});
  event_register ((Event){ 0, 1, init_0, {init_0_expr0}, ((int *)0), ((double *)0),
    "wall_heat.c", 101, "init"});
  event_register ((Event){ 0, 1, mov, {mov_expr0}, ((int *)0), ((double *)0),
    "wall_heat.c", 190, "mov"});
  event_register ((Event){ 0, 1, end, {end_expr0}, ((int *)0), ((double *)0),
    "wall_heat.c", 203, "end"});
  event_register ((Event){ 0, 1, set_dtmax, {set_dtmax_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 209, "set_dtmax"});
  event_register ((Event){ 0, 1, stability, {stability_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 211, "stability"});
  event_register ((Event){ 0, 1, vof, {vof_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 221, "vof"});
  event_register ((Event){ 0, 1, tracer_advection, {tracer_advection_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 222, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_advection_0, {tracer_advection_0_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/tracer.h", 42, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_diffusion, {tracer_diffusion_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 223, "tracer_diffusion"});
  event_register ((Event){ 0, 1, tracer_diffusion_0, {tracer_diffusion_0_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/tracer.h", 49, "tracer_diffusion"});
  event_register ((Event){ 0, 1, tracer_diffusion_1, {tracer_diffusion_1_expr0}, ((int *)0), ((double *)0),
    "wall_heat.c", 124, "tracer_diffusion"});
  event_register ((Event){ 0, 1, properties, {properties_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 230, "properties"});
  event_register ((Event){ 0, 1, properties_0, {properties_0_expr0}, ((int *)0), ((double *)0),
    "wall_heat.c", 51, "properties"});
  event_register ((Event){ 0, 1, advection_term, {advection_term_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 307, "advection_term"});
  event_register ((Event){ 0, 1, viscous_term, {viscous_term_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 337, "viscous_term"});
  event_register ((Event){ 0, 1, acceleration, {acceleration_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 373, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_0, {acceleration_0_expr0}, ((int *)0), ((double *)0),
    "wall_heat.c", 119, "acceleration"});
  event_register ((Event){ 0, 1, projection, {projection_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 416, "projection"});
  event_register ((Event){ 0, 1, end_timestep, {end_timestep_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 431, "end_timestep"});
  event_register ((Event){ 0, 1, adapt, {adapt_expr0}, ((int *)0), ((double *)0),
    "/home/dai/software/basilisk/src/navier-stokes/centered.h", 441, "adapt"});
  event_register ((Event){ 0, 1, adapt_0, {adapt_0_expr0}, ((int *)0), ((double *)0),
    "wall_heat.c", 199, "adapt"});
  _attribute = (_Attributes *) pcalloc (datasize/sizeof(double), sizeof (_Attributes), __func__, __FILE__, __LINE__);
  all = (scalar *) pmalloc (sizeof (scalar)*17,__func__, __FILE__, __LINE__);
  for (int i = 0; i < 16; i++)
    all[i].i = i;
  all[16].i = -1;
  set_fpe();
  quadtree_methods();
  init_face_vector ((vector){{14},{15}}, "muv");
  init_face_vector ((vector){{12},{13}}, "av");
  init_scalar ((scalar){11}, "b");
  init_face_vector ((vector){{9},{10}}, "uf");
  init_scalar ((scalar){8}, "pf");
  init_vector ((vector){{6},{7}}, "g");
  init_vector ((vector){{4},{5}}, "u");
  init_scalar ((scalar){3}, "p");
  embed = new_bid();
  init_face_vector ((vector){{1},{2}}, "fs");
  init_scalar ((scalar){0}, "cs");
  init_const_scalar ((scalar){_NVARMAX+5}, "zeroc",  0.);
  init_const_scalar ((scalar){_NVARMAX+4}, "unity",  1.);
  init_const_vector ((vector){{_NVARMAX+2},{_NVARMAX+3}}, "unityf", (double []) {1.,1.,1.});
  init_const_vector ((vector){{_NVARMAX+0},{_NVARMAX+1}}, "zerof", (double []) {0.,0.,0.});
  _set_boundary0();
  _set_boundary1();
  _set_boundary2();
  _set_boundary3();
  _set_boundary4();
  _set_boundary5();
  _set_boundary6();
  _set_boundary7();
  _set_boundary8();
  _set_boundary9();
  _set_boundary10();
  _set_boundary11();
}
