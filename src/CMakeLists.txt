

set(SRCS
    colormatrix.cpp
    )
    
set(HEADERS
    colormatrix.hpp
    )

    
add_library(essc_lib ${SRCS} ${HEADERS})
target_include_directories(essc_lib
    PUBLIC
    ./
)
target_link_libraries(essc_lib kmcapi sshash_static z)