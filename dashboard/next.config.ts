import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  // Fully static output — deploys to Vercel (or anywhere) with no server needed.
  output: "export",
};

export default nextConfig;
