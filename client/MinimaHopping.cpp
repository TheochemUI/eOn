#include "MinimaHopping.h"
#include <stdio.h>
#include "Log.h"

MinimaHopping::MinimaHopping(Parameters *params)
{
	parameters = params;
}

MinimaHopping::~MinimaHopping(void)
{

}

void MinimaHopping::run(void)
{
	log("HELLO FROM MH\n");
}
