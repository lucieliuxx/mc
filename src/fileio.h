#ifndef FILEIO_H
#define FILEIO_H

void write_data(long int nstep, SYSTEM* S);

void write_CoM(long int nstep, SYSTEM *S);

void write_unwrapped_trajectory(long int nstep, SYSTEM *S);

void write_final_config(long int nstep, SYSTEM *S);

void outfile_param_write(SYSTEM* S);

void infile_param_read(char *inputscript, SYSTEM *S);
#endif
