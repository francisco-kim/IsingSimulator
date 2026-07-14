using Microsoft.AspNetCore.Components.Web;
using Microsoft.AspNetCore.Components.WebAssembly.Hosting;

using IsingWeb;
using IsingWeb.Services;

var builder = WebAssemblyHostBuilder.CreateDefault(args);
builder.RootComponents.Add<App>("#app");
builder.RootComponents.Add<HeadOutlet>("head::after");

builder.Services.AddSingleton<SimulationRunner>();
builder.Services.AddSingleton<SweepRunner>();
builder.Services.AddSingleton<LatticeRenderer>();

await builder.Build().RunAsync();
