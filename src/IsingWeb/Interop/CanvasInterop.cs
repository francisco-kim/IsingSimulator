using System.Runtime.InteropServices.JavaScript;
using System.Runtime.Versioning;

namespace IsingWeb.Interop;

/// <summary>
///     Fast-path canvas rendering: the RGBA buffer is handed to JavaScript as a
///     view over WASM memory (no serialisation) and blitted with putImageData.
/// </summary>
[SupportedOSPlatform("browser")]
public static partial class CanvasInterop
{
    public const string ModuleName = "isingCanvas";

    private static Task? _moduleImport;

    public static Task EnsureModuleLoadedAsync() =>
        _moduleImport ??= JSHost.ImportAsync(ModuleName, "../js/isingCanvas.js");

    [JSImport("initCanvas", ModuleName)]
    public static partial void InitCanvas(string canvasId, int width, int height);

    [JSImport("renderFrame", ModuleName)]
    public static partial void RenderFrame(
        string canvasId,
        [JSMarshalAs<JSType.MemoryView>] Span<byte> rgba);

    [JSImport("initArrowCanvas", ModuleName)]
    public static partial void InitArrowCanvas(string canvasId, int size);

    [JSImport("clearArrows", ModuleName)]
    public static partial void ClearArrows(string canvasId);

    [JSImport("drawVortexArrows", ModuleName)]
    public static partial void DrawVortexArrows(
        string canvasId,
        [JSMarshalAs<JSType.MemoryView>] Span<double> arrows,
        int count,
        double cellPixels);

    [JSImport("downloadCanvasPng", ModuleName)]
    public static partial void DownloadCanvasPng(string canvasId, string filename);

    [JSImport("downloadText", ModuleName)]
    public static partial void DownloadText(string filename, string content);
}
