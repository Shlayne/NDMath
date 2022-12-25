namespace nd
{
	template<typename T, uint32_t N>
	inline bool Loader<T, N>::Load(const FilepathIn& crFilepath) noexcept
	{
		Reset();

		if (crFilepath.empty())
			return Fail(LoaderErrorCode_EmptyFilepath);

		if (std::error_code errorCode; !std::filesystem::is_regular_file(crFilepath, errorCode))
			return Fail(LoaderErrorCode_FileDoesNotExist);

		std::string meshFile;
		try
		{
			std::ifstream file(crFilepath);

			if (!file.is_open())
				return Fail(LoaderErrorCode_FileCouldNotOpen);

			std::stringstream fileStream;
			fileStream << file.rdbuf();
			file.close();

			if (fileStream.fail() || file.fail())
				return Fail(LoaderErrorCode_FileLoadingFailed);

			meshFile = fileStream.str();
		}
		catch (...)
		{
			return Fail(LoaderErrorCode_FileLoadingFailed);
		}

		// TODO: parse

		return Succeed();
	}

	template<typename T, uint32_t N>
	inline LoaderErrorCode Loader<T, N>::GetErrorCode() const noexcept
	{
		return m_ErrorCode;
	}

	template<typename T, uint32_t N>
	inline const char* Loader<T, N>::GetErrorMessage() const noexcept
	{
		switch (m_ErrorCode)
		{
			case LoaderErrorCode_EmptyFilepath:
				return "The provided filepath was empty.";
			case LoaderErrorCode_FileDoesNotExist:
				return "The provided filepath was not associated with an existing file.";
			case LoaderErrorCode_FileLoadingFailed:
				return "The provided filepath was associated with an existing file, but it failed to load.";
		}

		return nullptr;
	}

	template<typename T, uint32_t N>
	inline Loader<T, N>::MeshOut Loader<T, N>::GetOutputMesh() const noexcept
	{
		return m_MeshOut;
	}

	template<typename T, uint32_t N>
	inline void Loader<T, N>::Reset() noexcept
	{
		m_MeshOut.vertices.clear();
		m_MeshOut.indices.clear();
		m_ErrorCode = LoaderErrorCode_None;
	}

	template<typename T, uint32_t N>
	inline bool Loader<T, N>::Fail(LoaderErrorCode errorCode) noexcept
	{
		m_ErrorCode = errorCode;
		return false;
	}

	template<typename T, uint32_t N>
	inline bool Loader<T, N>::Succeed() noexcept
	{
		m_ErrorCode = LoaderErrorCode_None;
		return true;
	}
}
