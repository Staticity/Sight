#include <io/video_sources.hpp>
#include <utils/string_helpers.hpp>

#include <vector>

#include <mfapi.h>
#include <atlbase.h>
#include <strmif.h>
#include <propvarutil.h>
#include <dshow.h>

#pragma comment(lib, "strmiids.lib")

std::vector<sight::VideoSourceInfo> sight::EnumerateVideoSources()
{
    std::vector<VideoSourceInfo> sources;

    (void) CoInitializeEx(NULL, COINIT_APARTMENTTHREADED);
    HRESULT hr;
    ICreateDevEnum* pSysDevEnum = NULL;
    hr = CoCreateInstance(CLSID_SystemDeviceEnum, NULL, CLSCTX_INPROC_SERVER,
        IID_ICreateDevEnum, (void**)&pSysDevEnum);

    if (FAILED(hr))
    {
        return sources;
    }

    // Obtain a class enumerator for the video compressor category.
    IEnumMoniker* pEnumCat = NULL;
    hr = pSysDevEnum->CreateClassEnumerator(CLSID_VideoInputDeviceCategory, &pEnumCat, 0);

    if (hr == S_OK)
    {
        IBindCtx* pbc;
        hr = CreateBindCtx(NULL, &pbc);

        if (hr == S_OK)
        {
            // Enumerate the monikers.
            IMoniker* pMoniker = NULL;

            ULONG cFetched;
            while (pEnumCat->Next(1, &pMoniker, &cFetched) == S_OK)
            {
                IPropertyBag* pPropBag;
                hr = pMoniker->BindToStorage(pbc, 0, IID_IPropertyBag,
                    (void**)&pPropBag);

                if (SUCCEEDED(hr))
                {
                    VideoSourceInfo source;
                    source.index = sources.size();

                    // To retrieve the filter's friendly name, do the following:

                    std::vector<wchar_t*> fields = {
                        L"FriendlyName",
                        L"Description",
                        L"DevicePath"
                    };

                    for (int i = 0; i < fields.size(); ++i)
                    {
                        const auto& field = fields[i];

                        VARIANT varName;
                        VariantInit(&varName);
                        hr = pPropBag->Read(field, &varName, 0);

                        std::string value;
                        if (SUCCEEDED(hr))
                        {
                            value = duplicate_backslashes(string_cast(varName.bstrVal));
                        }
                        VariantClear(&varName);

                        switch (i)
                        {
                        case 0:
                            source.name = value;
                            break;
                        case 1:
                            source.description = value;
                            break;
                        case 2:
                            source.identifier = value;
                            break;
                        default:
                            break;
                        }
                    }

                    sources.push_back(source);

                    pPropBag->Release();
                }
                pMoniker->Release();
            }
            pEnumCat->Release();
        }

        pbc->Release();
    }
    pSysDevEnum->Release();
    return sources;
}
