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
	}

	includedirs {
		"src",
		"%{IncludeDir.gcem}",
		"%{IncludeDir.glm}",
	}

	filter "system:windows"
		systemversion "latest"
		defines "SYSTEM_WINDOWS"
		-- msvc doesn't provide __VA_OPT__ by default; this fixes that.
		usestdpreproc "On"
		enablemodules "Off"

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
