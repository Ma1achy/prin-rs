// APNG -> H.264 MP4, through ImageIO and AVFoundation.
//
// There is no ffmpeg on this machine and installing one to convert a file is a change to the
// user's system for a job macOS already does. ImageIO reads APNG frames; AVAssetWriter is the
// same encoder every other app on the machine uses.
//
// Usage: swift tools/apng2mp4.swift <in.png> <out.mp4> [fps] [mbps]
//
// The APNG stays the lossless record. This is lossy and inter-frame: a smooth ramp that bands in
// GIF's 256 colours survives here, but fine per-pixel structure is what H.264 spends bits on
// last, so read structure off the PNG and motion off this.

import AVFoundation
import CoreGraphics
import Foundation
import ImageIO

let args = CommandLine.arguments
guard args.count >= 3 else {
    FileHandle.standardError.write("usage: apng2mp4 <in.png> <out.mp4> [fps] [mbps]\n".data(using: .utf8)!)
    exit(2)
}
let inURL = URL(fileURLWithPath: args[1])
let outURL = URL(fileURLWithPath: args[2])
let fps = args.count > 3 ? Int32(args[3]) ?? 30 : 30
let mbps = args.count > 4 ? Double(args[4]) ?? 40.0 : 40.0

guard let src = CGImageSourceCreateWithURL(inURL as CFURL, nil) else {
    print("cannot open \(inURL.path)"); exit(1)
}
let n = CGImageSourceGetCount(src)
guard n > 1 else { print("only \(n) frame(s): ImageIO did not see an animation"); exit(1) }
guard let first = CGImageSourceCreateImageAtIndex(src, 0, nil) else { print("no frame 0"); exit(1) }
let (w, h) = (first.width, first.height)
print("in: \(n) frames, \(w)x\(h) -> \(fps) fps, \(mbps) Mb/s")

try? FileManager.default.removeItem(at: outURL)
let writer = try AVAssetWriter(outputURL: outURL, fileType: .mp4)
let settings: [String: Any] = [
    AVVideoCodecKey: AVVideoCodecType.h264,
    AVVideoWidthKey: w,
    AVVideoHeightKey: h,
    AVVideoCompressionPropertiesKey: [
        AVVideoAverageBitRateKey: Int(mbps * 1_000_000),
        AVVideoMaxKeyFrameIntervalKey: fps * 2,
    ],
]
let input = AVAssetWriterInput(mediaType: .video, outputSettings: settings)
input.expectsMediaDataInRealTime = false
let adaptor = AVAssetWriterInputPixelBufferAdaptor(
    assetWriterInput: input,
    sourcePixelBufferAttributes: [
        kCVPixelBufferPixelFormatTypeKey as String: Int(kCVPixelFormatType_32BGRA),
        kCVPixelBufferWidthKey as String: w,
        kCVPixelBufferHeightKey as String: h,
    ]
)
writer.add(input)
writer.startWriting()
writer.startSession(atSourceTime: .zero)

let cs = CGColorSpaceCreateDeviceRGB()
let started = Date()
for i in 0..<n {
    guard let img = CGImageSourceCreateImageAtIndex(src, i, nil) else { continue }
    while !input.isReadyForMoreMediaData { usleep(2000) }
    guard let pool = adaptor.pixelBufferPool else { print("no pixel buffer pool"); exit(1) }
    var pb: CVPixelBuffer?
    CVPixelBufferPoolCreatePixelBuffer(nil, pool, &pb)
    guard let buf = pb else { print("no pixel buffer"); exit(1) }
    CVPixelBufferLockBaseAddress(buf, [])
    if let ctx = CGContext(
        data: CVPixelBufferGetBaseAddress(buf),
        width: w, height: h, bitsPerComponent: 8,
        bytesPerRow: CVPixelBufferGetBytesPerRow(buf),
        space: cs,
        bitmapInfo: CGImageAlphaInfo.noneSkipFirst.rawValue | CGBitmapInfo.byteOrder32Little.rawValue
    ) {
        ctx.draw(img, in: CGRect(x: 0, y: 0, width: w, height: h))
    }
    CVPixelBufferUnlockBaseAddress(buf, [])
    adaptor.append(buf, withPresentationTime: CMTime(value: CMTimeValue(i), timescale: fps))
    if i % 100 == 0 || i == n - 1 {
        print("  frame \(i + 1)/\(n), \(String(format: "%.1f", -started.timeIntervalSinceNow))s")
    }
}
input.markAsFinished()
let done = DispatchSemaphore(value: 0)
writer.finishWriting { done.signal() }
done.wait()
if writer.status != .completed {
    print("writer failed: \(String(describing: writer.error))"); exit(1)
}
let bytes = (try? FileManager.default.attributesOfItem(atPath: outURL.path)[.size] as? Int) ?? 0
print("wrote \(outURL.path), \(String(format: "%.1f", Double(bytes ?? 0) / 1e6)) MB, \(n) frames")
