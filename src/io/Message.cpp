#include "Message.h"

#ifndef NDEBUG
Message::Verbosity Message::verbosity = QUIET;
tbb::mutex Message::mut;
#endif
