/**
 * \file CorrectScientificNotificationMain.h
 * 28/5/2010 LB Initial implementation
 * 
 */
#include "StringTools.h"

int main( int argc, const char* argv[] )
{
	if (argc < 2)
	{
		std::cout << "Usage: correctScientificNotation filename" << std::endl;
		return 1;
	}

	std::string filename = argv[1];
	correctScientificNotation(filename);
}