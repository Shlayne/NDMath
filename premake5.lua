include "Dependencies/premake/Custom/solutionitems.lua"

workspace "NDMath"
	architecture "x86_64"
	startproject "NDMath"

	configurations { "Debug", "Release" }

	solutionitems {
		-- Visual Studio
		".editorconfig",

		-- Git
		".gitignore",
		".gitattributes",

		-- Scripts
		"Scripts/GenerateProjects.bat",

		-- Lua Scripts
		"premake5.lua",
		"Dependencies/Dependencies.lua",
		"Dependencies/premake/Custom/solutionitems.lua",
		"Dependencies/premake/Custom/notyetimplemented.lua",

		-- Misc
		"README.md"
	}

	flags {
		"MultiProcessorCompile"
	}

OutputDir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"

include "Dependencies/premake/Custom/notyetimplemented.lua"
include "Dependencies/Dependencies.lua"

-- Add any projects here with 'include "__PROJECT_NAME__"'
include "NDMath"
