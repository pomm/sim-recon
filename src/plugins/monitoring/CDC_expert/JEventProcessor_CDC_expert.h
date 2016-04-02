// $Id$
//
//    File: JEventProcessor_CDC_expert.h
// Created: Wed Oct 22 2014
// Creator: Naomi Jarvis
//

#ifndef _JEventProcessor_CDC_expert_
#define _JEventProcessor_CDC_expert_

#include <JANA/JEventProcessor.h>


class JEventProcessor_CDC_expert:public jana::JEventProcessor{
 public:
  JEventProcessor_CDC_expert();
  ~JEventProcessor_CDC_expert();
  const char* className(void){return "JEventProcessor_CDC_expert";}

	bool initialized_histograms;

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_CDC_expert_

