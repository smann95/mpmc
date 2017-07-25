#pragma once
#ifndef OUTPUT_H
#define OUTPUT_H

class Output
{
public:
	Output();
	~Output();

	static void err( const char *msg);
	// Immediately displays a message to stdout

	static void out( const char *msg);
	// Immediately displays a message to stderr

	static char * make_filename( const char * basename, int fileno );
	// Forms a numbered filename, based on a common basename

	static double calctimediff(struct timeval a, struct timeval b);

	static int GetTimeOfDay( struct timeval *tv );
	
};



#endif // OUTPUT_H