#ifndef MINIMAHOPPING_H
#define MINIMAHOPPING_H
#include "Parameters.h"

class MinimaHopping
{
	public:
		MinimaHopping(Parameters *params);
		~MinimaHopping(void);
		void run(void);

	private:
		Parameters *parameters;
};

#endif
