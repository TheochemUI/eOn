#include "MinimaHoppingJob.h"
#include "MinimaHopping.h"
#include "Log.h"

MinimaHoppingJob::MinimaHoppingJob(Parameters *params)
{
	parameters = params;
}

MinimaHoppingJob::~MinimaHoppingJob(void)
{

}

std::vector<std::string> MinimaHoppingJob::run(void)
{
	log("hello %f\n", 3.14);
	MinimaHopping mh = MinimaHopping(parameters);
	mh.run();
    std::vector<std::string> returnFiles;
	return returnFiles;
}
