%% Progress bar compatible with parfor loops.
%
% queue = parforWaitbar(total, message) opens a waitbar and returns a
% parallel.pool.DataQueue.
%
% Workers call `send(queue, x)` once per completed iteration; the bar ticks on
% the client via an afterEach callback and auto-closes once `total` ticks have
% arrived.
%
% Based on https://vladislav-morozov.github.io/blog/simulations/tools/2024-11-11-simple-parfor-progress-bar/
%
function queue = parforWaitbar(total, message)

  h = waitbar(0, message);
  queue = parallel.pool.DataQueue;
  done = 0;
  afterEach(queue, @(~) tick());

  function tick()
    done = done + 1;
    waitbar(done / total, h, sprintf('%s (%d/%d)', message, done, total));
    if done >= total
      close(h);
    end
  end
end
