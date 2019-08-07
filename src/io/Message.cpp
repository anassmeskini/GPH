#include "Message.h"

Message::Verbosity Message::verbosity = RELEASE;

#ifndef NDEBUG
tbb::mutex Message::mut;
#endif
