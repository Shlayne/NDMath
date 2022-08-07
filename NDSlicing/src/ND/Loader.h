#pragma once

#include "ND/Matrix.h"
#include "ND/Mesh.h"
#include <filesystem>
#include <fstream>

namespace nd
{
	using LoaderErrorCode = uint8_t;
	enum : LoaderErrorCode
	{
		LoaderErrorCode_None = 0,
		LoaderErrorCode_EmptyFilepath = 1,
		LoaderErrorCode_FileDoesNotExist = 2,
		LoaderErrorCode_FileCouldNotOpen = 3,
		LoaderErrorCode_FileLoadingFailed = 4,
	};

	// Loads meshes from files.
	template<typename T, uint32_t N>
	class Loader
	{
	public:
		static_assert(N > 0, "Loader must have mesh with enough dimensions to load.");
	public:
		using FilepathIn = std::filesystem::path;
		using MeshOut = Mesh<T, N>;
	public:
		// Returns true if a mesh could be loaded, false if not or an error occurred.
		bool Load(const FilepathIn& crFilepath) noexcept;

		// Gets the error code of the most recent mesh loaded.
		LoaderErrorCode GetErrorCode() const noexcept;

		// Gets the error message of the most recent mesh loaded.
		const char* GetErrorMessage() const noexcept;

		// Gets the output mesh of the most recent mesh loaded.
		MeshOut GetOutputMesh() const noexcept;

		// Releases memory allocated by the most recent mesh loaded.
		// This is also done automatically at the beginning of Load, so
		// there's no need to call this in between successive Load calls.
		// But, after you're done using this Loader, calling this is recommended.
		void Reset() noexcept;
	private:
		bool Fail(LoaderErrorCode errorCode) noexcept;
		bool Succeed() noexcept;
	private:
		MeshOut m_MeshOut;
		LoaderErrorCode m_ErrorCode = LoaderErrorCode_None;
	};
}

#include "Loader.inl"
