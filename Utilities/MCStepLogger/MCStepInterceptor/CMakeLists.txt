# @author Sandro Wenzel
# @brief  cmake setup for module Utilities/MCStepLogger/MCStepLoggerImpl

set(MODULE_NAME "MCStepInterceptor")

# we do not depend on anything ... so using an empty bucket
set(MODULE_BUCKET_NAME mclogger_bucket)
O2_SETUP(NAME ${MODULE_NAME})

set(SRCS
    src/MCStepLogger.cxx
   )

set(LIBRARY_NAME ${MODULE_NAME})
set(BUCKET_NAME ${MODULE_BUCKET_NAME})

O2_GENERATE_LIBRARY()



