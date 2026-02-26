find_package(raylib CONFIG REQUIRED)

add_executable (FluidSimulation "FluidSimulation.cpp"
	"FluidSimulation.h" "main.cpp")

target_link_libraries(FluidSimulation PRIVATE raylib)

if (MSVC)
    target_compile_options(FluidSimulation PRIVATE /W4 /permissive-)
endif()
