project "NDMath"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++latest"
	cdialect "C17"
	staticruntime "On"

	targetdir ("%{wks.location}/bin/" .. OutputDir .. "/%{prj.name}")
	objdir ("%{wks.location}/bin-int/" .. OutputDir .. "/%{prj.name}")

	files {
		"src/**.h",
		"src/**.c",
		"src/**.hpp",
		"src/**.cpp",
		"src/**.inl",
		"src/**.ixx"
	}

	includedirs {
		"src",
		"%{IncludeDir.gcem}",
		"%{IncludeDir.glm}",
	}

	filter "system:windows"
		systemversion "latest"
		defines "SYSTEM_WINDOWS"

		-- Modules are OP.
		scanformoduledependencies "True"
		enablemodules "On"

		-- msvc doesn't provide __VA_OPT__ by default; this fixes that.
		usestdpreproc "On"

		-- These two are required because visual studio precompiled their module
		-- ifc's with dynamic linking and imprecise floating point operations.
		staticruntime "Off"
		floatingpoint "None"

	filter "configurations:Debug"
		runtime "Debug"
		optimize "Debug"
		symbols "Full"
		defines "CONFIG_DEBUG"

	filter "configurations:Release"
		runtime "Release"
		optimize "On"
		symbols "On"
		defines "CONFIG_RELEASE"
