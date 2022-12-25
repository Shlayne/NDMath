project "NDMath"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++latest"
	cdialect "C17"
	staticruntime "Off"
	floatingpoint "None"

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
		-- Add any project source directories here.
		"src",
		-- "%{wks.location}/__PROJECT_NAME__/src",

		-- Add any dependency includes here.
		"%{IncludeDir.gcem}",
		"%{IncludeDir.glm}",
	}

	-- Add any links dependency libs via their project names here.
	links {
		--	"__PROJECT_NAME__"
	}

	filter "system:windows"
		systemversion "latest"
		buildoptions "/wd5105" -- Until Microsoft updates Windows 10 to not have terrible code (aka never), this must be here to prevent a warning.
		defines "SYSTEM_WINDOWS"
		scanformoduledependencies "True"
		enablemodules "On"
		usestdpreproc "On"

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
