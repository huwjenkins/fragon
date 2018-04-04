cmake_minimum_required(VERSION 2.8)
project(fragon NONE)

install(
	DIRECTORY
		fragon
	DESTINATION
		lib/py2
)

install(
	DIRECTORY
		include
	DESTINATION
		share/fragon
)

install(
	PROGRAMS
		bin/fragon
	DESTINATION
		bin
)

if (WIN32)
    install(
        PROGRAMS
            bin/fragon.bat
        DESTINATION
            bin
    )
endif()